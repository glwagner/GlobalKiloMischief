# Astronomical tides in the stratified global ocean--sea ice model (include after setup.jl): an equilibrium tidal
# potential of M2, S2, K1 and O1, spherical-harmonic self-attraction and loading with PREM load Love numbers updated
# every time step, thickness-weighted energy-conserving Coriolis on partial cells, a biharmonic viscosity, quadratic
# bottom drag, an optional linear topographic wave drag, bathymetry regridded 11000 m deep, and an online harmonic fit
# of the sea surface height that survives restarts. The calibration behind these choices is in GlobalKiloMischief's
# glw/tides-1deg branch (tides/analysis/matrix/CALIBRATION.md).

using SpeedyTransforms, ConservativeRegridding  # load SphericalHarmonicSelfAttractionAndLoading
using Oceananigans.Coriolis: ThicknessWeightedEnergyConserving
using Oceananigans.Architectures: on_architecture
using Oceananigans.Operators: Δxᶜᶜᶜ, Δyᶜᶜᶜ
using JLD2: jldsave, load
using LinearAlgebra: inv
using KernelAbstractions: @kernel, @index
using Oceananigans.BuoyancyFormulations: ∂z_b
using Oceananigans.Grids: inactive_cell
using Oceananigans.Utils: launch!

const Adapt = NumericalEarth.Oceans.Adapt
const XDirection = NumericalEarth.Oceans.XDirection
const YDirection = NumericalEarth.Oceans.YDirection

#####
##### Equilibrium tidal potential
#####

"""
    EquilibriumTidalPotential(epoch; ramp = 2days, g = 9.80665, love_factor = 0.69)

The equilibrium tidal potential of M2, S2, K1 and O1 (Cartwright & Tayler amplitudes 0.242334, 0.112743, 0.141565 and
0.100661 m) times `love_factor` = 1 + k₂ - h₂, with astronomical arguments from the mean longitudes of the Moon and Sun
and nodal corrections f and u evaluated at the model time, `t` seconds after `epoch`, ramped up as (1 - cos(π t / ramp)) / 2
over `ramp`. Φ = -g η_eq, so that -∇Φ forces the flow towards the equilibrium tide
η_eq = Σ f A cos²φ cos(V + u + 2λ) (semidiurnal) + Σ f A sin 2φ cos(V + u + λ) (diurnal).
"""
struct EquilibriumTidalPotential{FT} <: Function
    centuries :: FT  # Julian centuries from J2000 to the epoch
    hours :: FT      # UT hour of the epoch
    ramp :: FT
    g :: FT
end

function EquilibriumTidalPotential(epoch; ramp = 2days, g = 9.80665, love_factor = 0.69)
    centuries = (Dates.datetime2julian(epoch) - 2451545.0) / 36525
    hours = Dates.hour(epoch) + Dates.minute(epoch) / 60
    return EquilibriumTidalPotential(centuries, hours, Float64(ramp), g * love_factor)
end

@inline function (Φ::EquilibriumTidalPotential)(λ, φ, t)
    T = Φ.centuries + t / (36525 * 86400)
    s = 218.3164477 + 481267.88123421 * T  # mean longitude of the Moon
    h = 280.4664567 + 36000.76983 * T      # mean longitude of the Sun
    N = 125.04452 - 1934.136261 * T        # longitude of the Moon's ascending node
    τ = 180 + 15 * (Φ.hours + t / 3600)    # hour angle of the mean Sun at Greenwich, plus 180ᵒ

    f_M2 = 1.0004 - 0.0373 * cosd(N) + 0.0002 * cosd(2N)
    u_M2 = -2.14 * sind(N)
    f_K1 = 1.0060 + 0.1150 * cosd(N) - 0.0088 * cosd(2N)
    u_K1 = -8.86 * sind(N) + 0.68 * sind(2N)
    f_O1 = 1.0089 + 0.1871 * cosd(N) - 0.0147 * cosd(2N)
    u_O1 = 10.80 * sind(N) - 1.34 * sind(2N) + 0.19 * sind(3N)

    semidiurnal = cosd(φ)^2
    diurnal = sind(2φ)
    η = 0.242334 * f_M2 * semidiurnal * cosd(2τ - 2s + 2h + u_M2 + 2λ) +
        0.112743 * semidiurnal * cosd(2τ + 2λ) +
        0.141565 * f_K1 * diurnal * cosd(τ + h - 90 + u_K1 + λ) +
        0.100661 * f_O1 * diurnal * cosd(τ - 2s + h + 90 + u_O1 + λ)

    ramp = ifelse(t < Φ.ramp, (1 - cospi(t / Φ.ramp)) / 2, one(t))
    return - Φ.g * ramp * η
end

#####
##### Bathymetry
#####

"""
    tidal_minimum_depth(Nz; depth = 11000, at_least = 5)

The depth of the shallowest face of the `Nz`-level vertical grid that is at least `at_least` meters deep, so that the
shallowest column holds whole surface cells.
"""
function tidal_minimum_depth(Nz; depth = 11000, at_least = 5)
    grid = RectilinearGrid(CPU(); size = Nz, z = ExponentialDiscretization(Nz, -depth, 0; scale = 1100),
                           topology = (Flat, Flat, Bounded))
    return -maximum(z for z in znodes(grid, Face()) if z ≤ -at_least)
end

"""
    tidal_bottom_height(Nx, Ny; minimum_depth)

ETOPO2022 bottom height regridded on an 11000 m deep grid, so that trenches are not clipped, with the one-cell-wide
channels of the Lincoln Sea (which resonate at M2) filled.
"""
function tidal_bottom_height(Nx, Ny; minimum_depth)
    zb = global_bottom_height(Nx, Ny; halo = (7, 7, 4), minimum_depth, regrid_depth = 11000)
    lincoln_sea = (λ = (280, 320), φ = (78, 84))
    filled = fill_one_cell_channels!(zb, lincoln_sea)
    @info "Filled $filled one-cell-wide channel cells in the Lincoln Sea on the $Nx × $Ny grid"
    return Array(interior(zb))
end

# Wet cells with dry neighbors on both sides in x or in y, inside the region, are made land until none are left.
function fill_one_cell_channels!(zb, region)
    grid = zb.grid
    Nx, Ny = size(grid, 1), size(grid, 2)
    z = view(interior(zb), :, :, 1)
    inside(i, j) = (λ = mod(λnode(i, j, 1, grid, Center(), Center(), Center()), 360);
                    φ = φnode(i, j, 1, grid, Center(), Center(), Center());
                    region.λ[1] ≤ λ ≤ region.λ[2] && region.φ[1] ≤ φ ≤ region.φ[2])
    dry(i, j) = z[mod1(i, Nx), clamp(j, 1, Ny)] ≥ 0
    filled = 0
    for pass in 1:20
        cells = [(i, j) for j in 2:Ny-1, i in 1:Nx if !dry(i, j) && inside(i, j) &&
                 ((dry(i-1, j) && dry(i+1, j)) || (dry(i, j-1) && dry(i, j+1)))]
        isempty(cells) && break
        for (i, j) in cells
            @info @sprintf("  one-cell channel at (%d, %d), λ = %.1f, φ = %.1f, %.0f m deep", i, j,
                           λnode(i, j, 1, grid, Center(), Center(), Center()), φnode(i, j, 1, grid, Center(), Center(), Center()), -z[i, j])
            z[i, j] = 0
        end
        filled += length(cells)
    end
    return filled
end

#####
##### Linear topographic wave drag
#####

"""
    WaveDrag(direction, rate)

Linear drag -r u (Jayne & St. Laurent 2001), with r at cell centers averaged to velocity faces.
"""
struct WaveDrag{D, R}
    direction :: D
    rate :: R
end

Adapt.@adapt_structure WaveDrag

@inline (d::WaveDrag{<:XDirection})(i, j, k, grid, clock, fields) = @inbounds - (d.rate[i-1, j, 1] + d.rate[i, j, 1]) / 2 * fields.u[i, j, k]
@inline (d::WaveDrag{<:YDirection})(i, j, k, grid, clock, fields) = @inbounds - (d.rate[i, j-1, 1] + d.rate[i, j, 1]) / 2 * fields.v[i, j, k]

"""
    wave_drag_rate(grid, bottom_height; χ, Δt)

r = χ C / H on the whole grid, with C = κ h² N_b √(1 - f²/ω²) / 2 on the 0.5ᵒ tiles of `tides/wavedrag_halfdeg_C.bin`
(zero where a tile is shallower than 100 m), replaced by r / (1 + r Δt) so that the forward step damps as backward Euler.
"""
function wave_drag_rate(grid, bottom_height; χ, Δt)
    C = reshape(reinterpret(Float64, read(joinpath(@__DIR__, "wavedrag_halfdeg_C.bin"))), 720, 312)  # 77.75ᵒS to 77.75ᵒN
    cpu_grid = TripolarGrid(CPU(); size = (size(bottom_height, 1), size(bottom_height, 2), 1), z = (-1, 0))
    rate = zeros(size(bottom_height, 1), size(bottom_height, 2))
    for j in axes(rate, 2), i in axes(rate, 1)
        H = -bottom_height[i, j, 1]
        λ = mod(λnode(i, j, 1, cpu_grid, Center(), Center(), Center()), 360)
        φ = φnode(i, j, 1, cpu_grid, Center(), Center(), Center())
        tile_i = clamp(floor(Int, λ / 0.5) + 1, 1, 720)
        tile_j = floor(Int, (φ + 78) / 0.5) + 1
        c = 1 ≤ tile_j ≤ 312 ? C[tile_i, tile_j] : 0.0
        r = H > 0 ? χ * c / H : 0.0
        rate[i, j] = r / (1 + r * Δt)
    end
    deep = rate[-bottom_height[:, :, 1] .> 1000]
    @info @sprintf("Wave drag χ = %g: median rate %.3f day⁻¹ where deeper than 1000 m, maximum %.3f day⁻¹", χ,
                   86400 * sort(deep)[end ÷ 2], 86400 * maximum(rate))
    r = Field{Center, Center, Nothing}(grid)
    set!(r, reshape(rate, size(rate)..., 1))
    fill_halo_regions!(r)
    return r
end

"""
    subgrid_roughness_over_2H(grid, bottom_height, cells_per_degree)

S / 2H on the whole grid, with S the roughness of topography finer than twice the grid spacing from
`tides/subgrid_roughness.jld2` (written by `tides/subgrid_wave_drag.jl`), on 1ᵒ tiles.
"""
function subgrid_roughness_over_2H(grid, bottom_height, cells_per_degree)
    S = load(joinpath(@__DIR__, "subgrid_roughness.jld2"), cells_per_degree == 6 ? "S_sixth" : "S_twelfth")
    cpu_grid = TripolarGrid(CPU(); size = (size(bottom_height, 1), size(bottom_height, 2), 1), z = (-1, 0))
    S_over_2H = zeros(size(bottom_height, 1), size(bottom_height, 2))
    for j in axes(S_over_2H, 2), i in axes(S_over_2H, 1)
        H = -bottom_height[i, j, 1]
        λ = λnode(i, j, 1, cpu_grid, Center(), Center(), Center())
        φ = φnode(i, j, 1, cpu_grid, Center(), Center(), Center())
        tile_i = clamp(floor(Int, mod(λ + 180, 360)) + 1, 1, 360)
        tile_j = clamp(floor(Int, φ + 90) + 1, 1, 180)
        S_over_2H[i, j] = H > 0 ? S[tile_i, tile_j] / 2H : 0.0
    end
    field = Field{Center, Center, Nothing}(grid)
    set!(field, reshape(S_over_2H, size(S_over_2H)..., 1))
    return field
end

"""
    set_scale_aware_wave_drag!(r, ocean_model, S_over_2H; χ = 1, Δt, ω = 2π / 12.4206012hours, layers = 3)

Fill the wave drag rate r = χ S N_b √(1 - f²/ω²) / 2H (Johnston, Shakespeare & Constantinou 2025, Eq. 4, restricted to
unresolved topography), with N_b the buoyancy frequency averaged over the `layers` lowest interfaces of each column
of the current state, zero where f > ω (no propagating waves), and r replaced by r / (1 + r Δt).
"""
function set_scale_aware_wave_drag!(r, ocean_model, S_over_2H; χ = 1, Δt, ω = 2π / 12.4206012hours, layers = 3)
    grid = ocean_model.grid
    launch!(architecture(grid), grid, :xy, _scale_aware_wave_drag!, r, grid, S_over_2H, ocean_model.buoyancy,
            ocean_model.tracers, χ, Δt, ω, layers)
    fill_halo_regions!(r)
    rate = filter(>(0), Array(interior(r)))
    isempty(rate) || @info @sprintf("Scale-aware wave drag χ = %g: median rate %.4f day⁻¹, maximum %.3f day⁻¹, %d cells", χ,
                                    86400 * sort(rate)[end ÷ 2], 86400 * maximum(rate), length(rate))
    return nothing
end

@kernel function _scale_aware_wave_drag!(r, grid, S_over_2H, buoyancy, tracers, χ, Δt, ω, layers)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)
    bottom = Nz + 1
    for k in Nz:-1:1
        bottom = ifelse(inactive_cell(i, j, k, grid), bottom, k)
    end
    N² = zero(grid)
    n = 0
    for k in bottom+1:min(bottom + layers, Nz)
        N² += ∂z_b(i, j, k, grid, buoyancy, tracers)
        n += 1
    end
    N = sqrt(max(0, N² / max(n, 1)))
    f = 2 * 7.292115e-5 * sind(φnode(i, j, 1, grid, Center(), Center(), Center()))
    rate = χ * @inbounds(S_over_2H[i, j, 1]) * N * sqrt(max(0, 1 - f^2 / ω^2))
    @inbounds r[i, j, 1] = rate / (1 + rate * Δt)
end

#####
##### Ocean configuration
#####

@inline Δ²ᶜᶜᶜ(i, j, k, grid) = (Δx = Δxᶜᶜᶜ(i, j, k, grid); Δy = Δyᶜᶜᶜ(i, j, k, grid); 2Δx^2 * Δy^2 / (Δx^2 + Δy^2))
# Evaluated where needed (from the cell-center spacings) rather than stored, which would take a 3D field.
@inline biharmonic_A(i, j, k, grid, ℓx, ℓy, ℓz, clock, fields, τ) = Δ²ᶜᶜᶜ(i, j, k, grid)^2 / τ

"""
    tidal_ocean_kw(grid, bottom_height; max_Δt, wave_drag = :none, wave_drag_coefficient = 1,
                   cells_per_degree, biharmonic = false)

`ocean_simulation` keywords for tides: the equilibrium potential of M2, S2, K1 and O1 from `start_date`, ramped over 2
days; spherical-harmonic SAL (degree 40, PREM load Love numbers); `ThicknessWeightedEnergyConserving` Coriolis; CATKE,
and if `biharmonic` a biharmonic viscosity A = Δ⁴/τ, τ = max(1 day, 64 max_Δt) (the explicit stability bound), which
Oceananigans#6451 makes dissipative on partial cells; quadratic bottom drag 0.003; and linear wave drag with coefficient
χ = `wave_drag_coefficient`: `wave_drag = :jsl` (Jayne & St. Laurent with the 0.5ᵒ C field) or `:scale_aware` (unresolved
roughness only; its rate is zero until `set_scale_aware_wave_drag!` fills it from the initial stratification).
"""
function tidal_ocean_kw(grid, bottom_height; max_Δt, wave_drag = :none, wave_drag_coefficient = 1,
                        cells_per_degree, biharmonic = false)
    closure = NumericalEarth.Oceans.default_ocean_closure()
    if biharmonic
        τ = max(1day, 64max_Δt)
        @info @sprintf("Biharmonic viscosity A = Δ⁴/τ, τ = %.2f days", τ / day)
        closure = (closure, HorizontalScalarBiharmonicDiffusivity(ν = biharmonic_A, discrete_form = true, parameters = τ))
    end

    love_numbers = NumericalEarth.Oceans.prem_load_love_numbers()
    self_attraction_and_loading = SphericalHarmonicSelfAttractionAndLoading(grid; truncation = 40, nlat_half = 24, love_numbers)

    forcing = if wave_drag == :jsl
        r = wave_drag_rate(grid, bottom_height; χ = wave_drag_coefficient, Δt = max_Δt)
        (u = WaveDrag(XDirection(), r), v = WaveDrag(YDirection(), r))
    elseif wave_drag == :scale_aware
        r = Field{Center, Center, Nothing}(grid)
        (u = WaveDrag(XDirection(), r), v = WaveDrag(YDirection(), r))
    else
        NamedTuple()
    end

    return (; barotropic_potential = EquilibriumTidalPotential(start_date), self_attraction_and_loading,
              coriolis = HydrostaticSphericalCoriolis(scheme = ThicknessWeightedEnergyConserving()),
              closure, bottom_drag_coefficient = 0.003, forcing)
end

"""
    update_self_attraction_and_loading_every_step!(ocean)

Replace the SAL update at every stage of the time stepper by one update per time step.
"""
function update_self_attraction_and_loading_every_step!(ocean)
    sal = ocean.callbacks[:self_attraction_and_loading].parameters
    delete!(ocean.callbacks, :self_attraction_and_loading)
    add_callback!(ocean, sim -> NumericalEarth.Oceans.compute_self_attraction_and_loading!(sim.model, sal),
                  IterationInterval(1); name = :self_attraction_and_loading_per_step)
    return nothing
end

#####
##### Online harmonic fit
#####

const fit_periods = (M2 = 12.4206012, S2 = 12.0, N2 = 12.65834751, K1 = 23.9344696, O1 = 25.8193417)  # hours

"""
    add_harmonic_fit!(simulation, η, directory; start = 20days, save_interval = 5days, first_window = 50days,
                      minimum_span = 15days)

Accumulate a least-squares fit of a mean and M2, S2, N2, K1 and O1 to hourly `η` from `start` on. Every `save_interval`
(with the checkpoints) each rank writes `harmonics_rank<r>.jld2` in `directory`: the sums, the sample times and, once
the samples span `minimum_span` (M2, S2 and N2 separate over about 15 days), the fitted complex amplitudes Z
(η ≈ Re Z exp(iωt), t in seconds from `start_date`), and once, at `first_window`, a copy
`harmonics_days<start>-<first_window>_rank<r>.jld2`. A restart reloads the sums and skips times already summed.
"""
function add_harmonic_fit!(simulation, η, directory; start = 20days, save_interval = 5days, first_window = 50days,
                           minimum_span = 15days)
    rank = MPI.Initialized() ? MPI.Comm_rank(MPI.COMM_WORLD) : 0
    path = joinpath(directory, "harmonics_rank$rank.jld2")
    ω = [2π / (3600T) for T in values(fit_periods)]
    sums = [Field{Center, Center, Nothing}(η.grid) for _ in 1:(2length(ω) + 1)]
    times = Float64[]
    if isfile(path)
        saved = load(path)
        append!(times, saved["times"])
        for (S, s) in zip(sums, saved["sums"])
            interior(S) .= on_architecture(architecture(S), s)
        end
        @info "Reloaded the harmonic fit from $path: $(length(times)) samples"
    end

    function accumulate(sim)
        t = time(sim)
        (t < start - 1 || (!isempty(times) && t ≤ last(times) + 1)) && return nothing
        push!(times, t)
        for (S, b) in zip(sums, [1.0; cos.(ω .* t); sin.(ω .* t)])
            interior(S) .+= b .* interior(η)
        end
        return nothing
    end

    function save(sim)
        isempty(times) && return nothing
        S = [Array(interior(F)) for F in sums]
        fit = (;)
        if last(times) - first(times) ≥ minimum_span
            B = hcat(ones(length(times)), [cos(w * t) for t in times, w in ω], [sin(w * t) for t in times, w in ω])
            G⁻¹ = inv(B' * B)
            c(n) = sum(G⁻¹[n, m] .* S[m] for m in eachindex(S))
            K = length(ω)
            Z = Dict(string(name) => ComplexF32.(c(1 + k) .- im .* c(1 + K + k)) for (k, name) in enumerate(keys(fit_periods)))
            fit = (; Z, mean = Float32.(c(1)))
        end
        temporary = path * ".partial"
        jldsave(temporary; sums = S, times, first_day = first(times) / day, last_day = last(times) / day, fit...)
        mv(temporary, path; force = true)
        window = joinpath(directory, @sprintf("harmonics_days%d-%d_rank%d.jld2", round(Int, start / day), round(Int, first_window / day), rank))
        time(sim) ≥ first_window && !isfile(window) && cp(path, window)
        return nothing
    end

    add_callback!(simulation, accumulate, TimeInterval(1hour); name = :harmonic_fit)
    add_callback!(simulation, save, TimeInterval(save_interval); name = :save_harmonic_fit)
    return nothing
end
