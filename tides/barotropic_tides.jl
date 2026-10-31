# Barotropic tides in a single-layer global ocean of uniform density forced by the astronomical tidal potential,
# for comparing treatments of self-attraction and loading (SAL) with each other and with altimetric atlases.
# The grid, bathymetry, bottom drag and biharmonic viscosity are as in internal_tides.jl. There is no wave drag,
# so the tides are expected to be too large. Output is the sea surface height every hour.
#
#   SELF_ATTRACTION_AND_LOADING=<none | scalar | spherical_harmonic | convolution> \
#       julia --project barotropic_tides.jl <CPU | CUDA> <cells per degree> <Δt in minutes> <days>
#
# The convolution (Chen, Wang, Arbic & Krasny 2026) is summed directly in CONVOLUTION_PRECISION (default Float32),
# with PREM load Love numbers and the self-cell term (CONVOLUTION_KERNEL=exact, the default), or as in the paper's MOM6
# code (CONVOLUTION_KERNEL=paper: asymptotic Love numbers, no self-cell term, kernel divided by 3).

using NumericalEarth
using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Tidejinks
using Dates
using JLD2
using Printf
const LinearAlgebra = Base.require(Base.PkgId(Base.UUID("37e2e46d-f89d-539d-b4ee-838fcccc9c8e"), "LinearAlgebra"))
import SPICE

arch_name = ARGS[1]
cells_per_degree = parse(Float64, ARGS[2])
Δt = parse(Float64, ARGS[3]) * minutes
stop_days = parse(Float64, ARGS[4])

arch_name == "CUDA" && @eval using CUDA
arch = arch_name == "CUDA" ? GPU() : CPU()

start_date = DateTime(2015, 1, 1)
end_date = start_date + Millisecond(round(Int, stop_days * 86400e3))

cd(@__DIR__)
mkpath("output")
Tidejinks.wrangle_spice_kernels("kernels.txt", "spice")
SPICE.furnsh("kernels.txt")

Nx, Ny = round(Int, 360 * cells_per_degree), round(Int, 170 * cells_per_degree)
depth = 6000
halo = (7, 7, 4)

cpu_grid = TripolarGrid(CPU(); size = (Nx, Ny, 1), halo, z = (-depth, 0))
bottom_height = regrid_bathymetry(cpu_grid; minimum_depth = 20, interpolation_passes = 3, major_basins = 1)

# With true depths (see MINIMUM_FRACTIONAL_CELL_HEIGHT below), isolated shallow cells (atolls and seamount tops, e.g. a
# 20 m cell at 172ᵒE 6ᵒN) surrounded by abyssal water go unstable. Cells whose wet neighbours are all more than
# ISOLATED_DEPTH_RATIO (default 5, 0 to disable) times deeper are made land, repeatedly until none are left.
isolated_ratio = parse(Float64, get(ENV, "ISOLATED_DEPTH_RATIO", "5"))
if isolated_ratio > 0
    let zb = interior(bottom_height, :, :, 1), removed = 0
        for pass in 1:10
            H = -Array(zb)
            n = 0
            for j in 2:Ny-1, i in 1:Nx
                H[i, j] > 0 || continue
                neighbours = (H[mod1(i-1, Nx), j], H[mod1(i+1, Nx), j], H[i, j-1], H[i, j+1])
                wet = filter(>(0), neighbours)
                if !isempty(wet) && all(>(isolated_ratio * H[i, j]), wet)
                    zb[i, j] = 0; n += 1
                end
            end
            removed += n
            n == 0 && break
        end
        @info "Made $removed isolated shallow cells land (ISOLATED_DEPTH_RATIO = $isolated_ratio)"
    end
end

underlying_grid = TripolarGrid(arch; size = (Nx, Ny, 1), halo, z = (-depth, 0))
zb = Field{Center, Center, Nothing}(underlying_grid)
set!(zb, Array(interior(bottom_height)))
fill_halo_regions!(zb)
# The single layer is 6000 m thick, so PartialCellBottom's default minimum_fractional_cell_height = 0.2 would make every
# column shallower than 1200 m 1200 m deep (wave speed, transport and bottom drag/Δz all wrong on shelves and slopes).
# MINIMUM_FRACTIONAL_CELL_HEIGHT (default 1e-3, i.e. 6 m < minimum_depth = 20 m) keeps the true depth everywhere.
ϵ = parse(Float64, get(ENV, "MINIMUM_FRACTIONAL_CELL_HEIGHT", "1e-3"))
grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(zb; minimum_fractional_cell_height = ϵ); active_cells_map = true)

Φ = Tidejinks.TidalPotential(grid, (start_date, end_date); reference_frame = "IAU_EARTH")
g = 9.80665

τ = parse(Float64, get(ENV, "BIHARMONIC_TIMESCALE", "0")) * days
@inline νhb(i, j, k, grid, τ) = Oceananigans.Operators.Azᶜᶜᶜ(i, j, k, grid)^2 / τ
ν = Field{Center, Center, Center}(grid)
τ > 0 && set!(ν, KernelFunctionOperation{Center, Center, Center}(νhb, grid, τ))
# BIHARMONIC_TIMESCALE (days; default 0 = off) and LAPLACIAN_VISCOSITY (m² s⁻¹, default 1000). With true depths, the biharmonic
# operator at steep partial-cell steps drives steady grid-scale η of tens to hundreds of metres, growing with ν.
νₗ = parse(Float64, get(ENV, "LAPLACIAN_VISCOSITY", "1000"))
closure = Tuple(c for c in (τ > 0 ? HorizontalScalarBiharmonicDiffusivity(; ν) : nothing,
                            νₗ > 0 ? HorizontalScalarDiffusivity(ν = νₗ) : nothing) if !isnothing(c))
closure = isempty(closure) ? nothing : length(closure) == 1 ? only(closure) : closure
free_surface = NumericalEarth.Oceans.default_free_surface(grid; fixed_Δt = Δt)

loading = get(ENV, "SELF_ATTRACTION_AND_LOADING", "none")
loading_suffix = loading
self_attraction_and_loading = if loading == "scalar"
    ScalarSelfAttractionAndLoading(parse(Float64, get(ENV, "SCALAR_LOADING_COEFFICIENT", "0.1")))
elseif loading == "spherical_harmonic"
    @eval using SpeedyTransforms, ConservativeRegridding
    truncation = parse(Int, get(ENV, "LOADING_TRUNCATION", "40"))
    SphericalHarmonicSelfAttractionAndLoading(grid; truncation, nlat_half = max(24, cld(truncation + 1, 2) + 1))
elseif loading == "convolution"
    FT = get(ENV, "CONVOLUTION_PRECISION", "Float32") == "Float64" ? Float64 : Float32
    if get(ENV, "CONVOLUTION_KERNEL", "exact") == "paper"
        loading_suffix = "convolution_paper"
        ConvolutionSelfAttractionAndLoading(grid; FT, love_numbers = AsymptoticLoadLoveNumbers(), self_term = false,
                                            earth_density = 3 * 5517)
    else
        loading_suffix = "convolution_prem"
        ConvolutionSelfAttractionAndLoading(grid; FT)
    end
else
    nothing
end
@info "Self-attraction and loading: $(isnothing(self_attraction_and_loading) ? "none" : summary(self_attraction_and_loading))"

# Linear topographic wave drag -χ C / H u (Jayne & St. Laurent 2001) with WAVE_DRAG_COEFFICIENT = χ > 0, where
# C = κ h² N_b √(1 - f²/ω²) / 2 on 0.5ᵒ tiles comes from egbert_ray/solve/scripts/prep_wavedrag.py (zero where the
# tile is shallower than 100 m), interpolated to cell centres. The rate χ C / H is capped at 1 / 2Δt for stability.
const Adapt = NumericalEarth.Oceans.Adapt
struct WaveDrag{D, R}
    direction :: D
    rate :: R
end
Adapt.@adapt_structure WaveDrag
@inline (d::WaveDrag{<:NumericalEarth.Oceans.XDirection})(i, j, k, grid, clock, fields) = @inbounds - (d.rate[i-1, j, 1] + d.rate[i, j, 1]) / 2 * fields.u[i, j, k]
@inline (d::WaveDrag{<:NumericalEarth.Oceans.YDirection})(i, j, k, grid, clock, fields) = @inbounds - (d.rate[i, j-1, 1] + d.rate[i, j, 1]) / 2 * fields.v[i, j, k]

χ = parse(Float64, get(ENV, "WAVE_DRAG_COEFFICIENT", "0"))
forcing = if χ > 0
    C = reshape(reinterpret(Float64, read("wavedrag_halfdeg_C.bin")), 720, 312) # 0.5ᵒ tiles, lat from -77.75 to 77.75
    function tile_value(λ, φ)
        i = clamp(floor(Int, mod(λ, 360) / 0.5) + 1, 1, 720)
        j = floor(Int, (φ + 78) / 0.5) + 1
        return 1 ≤ j ≤ 312 ? C[i, j] : 0.0
    end
    rate = zeros(Nx, Ny)
    for j in 1:Ny, i in 1:Nx
        H = -bottom_height[i, j, 1]
        λ = Oceananigans.Grids.λnode(i, j, 1, cpu_grid, Center(), Center(), Center())
        φ = Oceananigans.Grids.φnode(i, j, 1, cpu_grid, Center(), Center(), Center())
        rate[i, j] = H > 0 ? min(χ * tile_value(λ, φ) / H, 1 / 2Δt) : 0.0
    end
    r = Field{Center, Center, Nothing}(underlying_grid)
    set!(r, reshape(rate, Nx, Ny, 1))
    fill_halo_regions!(r)
    @info @sprintf("Wave drag χ = %g: rate χ C / H median %.3f day⁻¹ over cells deeper than 1000 m, capped in %d cells", χ,
                   86400 * sort(rate[-interior(bottom_height, :, :, 1) .> 1000])[end ÷ 2],
                   count(rate .== 1 / 2Δt))
    (u = WaveDrag(NumericalEarth.Oceans.XDirection(), r), v = WaveDrag(NumericalEarth.Oceans.YDirection(), r))
else
    NamedTuple()
end
wave_drag_suffix = χ > 0 ? @sprintf("_wavedrag%g", χ) : ""

# CORIOLIS_SCHEME = EnstrophyConserving (default), EnergyConserving, ActiveWeightedEnstrophyConserving,
# ActiveWeightedEnergyConserving or TriadScheme
scheme_name = get(ENV, "CORIOLIS_SCHEME", "EnstrophyConserving")
coriolis_scheme = getproperty(scheme_name ∈ ("EnstrophyConserving", "EnergyConserving") ? Oceananigans.Advection : Oceananigans.Coriolis,
                              Symbol(scheme_name))()
coriolis = HydrostaticSphericalCoriolis(scheme = coriolis_scheme)
coriolis_suffix = scheme_name == "EnstrophyConserving" ? "" : "_" * scheme_name

bottom_drag_coefficient = parse(Float64, get(ENV, "BOTTOM_DRAG_COEFFICIENT", "0.003"))
# MOMENTUM_ADVECTION = none (default: linear momentum, as WENO vector-invariant advection blows up next to 20 m cells) or default
momentum_advection = get(ENV, "MOMENTUM_ADVECTION", "none") == "none" ? (; momentum_advection = nothing) : (;)
ocean = ocean_simulation(grid; momentum_advection..., Δt, free_surface, barotropic_potential = Φ, closure, bottom_drag_coefficient, forcing, coriolis,
                         self_attraction_and_loading, equation_of_state = LinearEquationOfState(),
                         radiative_forcing = nothing, gravitational_acceleration = g, stop_time = stop_days * days)
set!(ocean.model, T = 10, S = 35)

η = ocean.model.free_surface.displacement
# RUN_TAG is appended to the output name.
name = @sprintf("barotropic_tides_%gdeg_dt%gmin_%s%s%s", 1 / cells_per_degree, Δt / minutes, loading_suffix,
                wave_drag_suffix * coriolis_suffix, get(ENV, "RUN_TAG", ""))
# FIT_START_DAY = d > 0: instead of hourly η, accumulate a least-squares fit of a mean and M2, N2, S2, K1, O1 to hourly
# η from day d to the end (N2 separates from M2 only over ≥ 27.6 days), and save the complex amplitudes
# Z (η ≈ Re Z exp(iωt), t from start_date) in output/<name>_harmonics.jld2, with η every 5 days for checking.
fit_start = parse(Float64, get(ENV, "FIT_START_DAY", "0")) * days
fit_periods = (M2 = 12.4206012, N2 = 12.65834751, S2 = 12.0, K1 = 23.9344696, O1 = 25.8193417) # hours
fit_ω = [2π / (3600T) for T in values(fit_periods)]
fit_times = Float64[]
fit_sums = [Field{Center, Center, Nothing}(grid) for _ in 1:(2length(fit_ω) + 1)]
η_snapshots = Dict{String, Array{Float32}}()
function accumulate_fit(sim)
    t = time(sim)
    t < fit_start - 1 && return nothing
    push!(fit_times, t)
    basis = [1.0; cos.(fit_ω .* t); sin.(fit_ω .* t)]
    for (S, b) in zip(fit_sums, basis)
        interior(S) .+= b .* interior(η)
    end
    return nothing
end
function save_fit(sim)
    t = fit_times
    A = hcat(ones(length(t)), [cos(w * τ) for τ in t, w in fit_ω], [sin(w * τ) for τ in t, w in fit_ω])
    G⁻¹ = LinearAlgebra.inv(A' * A)
    S = [Array(interior(F, :, :, 1)) for F in fit_sums]
    coefficient(n) = sum(G⁻¹[n, m] .* S[m] for m in eachindex(S))
    K = length(fit_ω)
    Z = Dict(string(c) => coefficient(1 + k) .- im .* coefficient(1 + K + k) for (k, c) in enumerate(keys(fit_periods)))
    jldsave(joinpath("output", name * "_harmonics.jld2"); Z, mean = coefficient(1), first_day = t[1] / 86400,
            last_day = t[end] / 86400, samples = length(t), η_snapshots)
    @info @sprintf("Saved %s_harmonics.jld2: fit of %d hourly samples, days %.2f–%.2f", name, length(t), t[1] / 86400, t[end] / 86400)
end
if fit_start > 0
    add_callback!(ocean, accumulate_fit, TimeInterval(1hours))
    add_callback!(ocean, sim -> (η_snapshots[@sprintf("day%03d", round(Int, time(sim) / 86400))] = Array(interior(η, :, :, 1))),
                  TimeInterval(5days))
else
    ocean.output_writers[:surface] = JLD2Writer(ocean.model, (; η); filename = joinpath("output", name),
                                                schedule = TimeInterval(1hours), overwrite_files = true)
end

wall_time = Ref(time_ns())
function progress(sim)
    elapsed = (time_ns() - wall_time[]) / 1e9
    ηc = Array(interior(η, :, :, 1))
    m = findmax(x -> isnan(x) ? Inf : abs(x), ηc)
    i, j = Tuple(m[2])
    @info @sprintf("%s: day %.2f, max|η| = %.2f m at (%d, %d) λ = %.1f φ = %.1f H = %.0f m, max|u| = %.2f m/s, wall %.1f s",
                   name, time(sim) / 86400, m[1], i, j, Oceananigans.Grids.λnode(i, j, 1, cpu_grid, Center(), Center(), Center()),
                   Oceananigans.Grids.φnode(i, j, 1, cpu_grid, Center(), Center(), Center()), -bottom_height[i, j, 1],
                   maximum(abs, sim.model.velocities.u), elapsed)
    flush(stderr)
    wall_time[] = time_ns()
end
add_callback!(ocean, progress, TimeInterval(parse(Float64, get(ENV, "PROGRESS_INTERVAL_DAYS", "1")) * days))

run!(ocean)
fit_start > 0 && save_fit(ocean)
