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
using Oceananigans.Operators: Δxᶜᶜᶜ, Δyᶜᶜᶜ, δxᶜᶜᶜ, δyᶜᶜᶜ, δxᶠᶠᶜ, δyᶠᶠᶜ, Δxᶠᶠᶜ, Δyᶠᶠᶜ, ℑxyᶜᶜᵃ
using Oceananigans.Operators: Axᶠᶜᶜ, Ayᶜᶠᶜ, Vᶠᶜᶜ, Vᶜᶠᶜ, Δzᶜᶜᶜ, ℑxᶠᵃᵃ, ℑyᵃᶠᵃ, ℑxᶜᵃᵃ, ℑyᵃᶜᵃ
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

# MINIMUM_DEPTH (m, default 20). EXTERNAL_DEPTH = mom6: Chen et al.'s MOM6 p04arc9 depth (aw, moreocean, ocnpct0p5) and JSL_RH,
# regridded onto this grid by analysis/diag/regrid_mom6.jl (mom6_p04arc9_on_<Nx>.jld2), instead of regrid_bathymetry. Wet
# cells are deepened to MINIMUM_DEPTH (MOM6 itself keeps cells down to 5 m above sea level wet, with thickness D + η), and the
# single layer is made as deep as the deepest cell (MOM6 has 9314 m). Without EXTERNAL_DEPTH the bathymetry is regridded onto a
# 11000 m deep grid first and the layer is then made as deep as the deepest regridded cell, so PartialCellBottom (which clamps the
# bottom to the grid) no longer clips trenches and abyssal plains at the old 6000 m layer depth.
minimum_depth = parse(Float64, get(ENV, "MINIMUM_DEPTH", "20"))
external_depth = get(ENV, "EXTERNAL_DEPTH", "none")
external = external_depth == "mom6" ? load(@sprintf("mom6_p04arc9_on_%d.jld2", Nx)) : nothing
isnothing(external) || (depth = max(depth, ceil(maximum(external["depth"]))))
isnothing(external) && (depth = 11000)
cpu_grid = TripolarGrid(CPU(); size = (Nx, Ny, 1), halo, z = (-depth, 0))
bottom_height = if isnothing(external)
    regrid_bathymetry(cpu_grid; minimum_depth, interpolation_passes = 3, major_basins = 1)
else
    let b = Field{Center, Center, Nothing}(cpu_grid), D = external["depth"]
        set!(b, reshape(map(x -> x > 0 ? -max(x, minimum_depth) : 0.0, D), Nx, Ny, 1))
        fill_halo_regions!(b)
        @info @sprintf("External depth %s: %d wet cells, %d deepened to %g m, layer %g m deep", external_depth, count(D .> 0),
                       count(0 .< D .< minimum_depth), minimum_depth, depth)
        b
    end
end
isnothing(external) && (depth = ceil(-minimum(interior(bottom_height))))

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
# (The default is lowered to half of MINIMUM_DEPTH / layer depth if that is smaller.)
ϵ = parse(Float64, get(ENV, "MINIMUM_FRACTIONAL_CELL_HEIGHT", string(min(1e-3, minimum_depth / 2depth))))
grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(zb; minimum_fractional_cell_height = ϵ); active_cells_map = true)

g = parse(Float64, get(ENV, "GRAVITATIONAL_ACCELERATION", "9.80665")) # MOM6 (Chen et al.): 9.8
# M2_ONLY=true: equilibrium M2 alone (0.69 × 0.2423 m, nodal f and u fixed at 2015-01-16 as the scoring assumes for day 30-60 windows),
# η_eq = f A cos²φ cos(ωt + V₀ + u + 2λ), and the online fit is of M2 and a mean only, as in the Chen et al. setup.
m2_only = get(ENV, "M2_ONLY", "false") == "true"
const M2_G = g; const M2_ω = 2π / (12.4206012 * 3600)
const M2_A = let jc = (Dates.datetime2julian(DateTime(2015, 2, 15)) - 2451545.0) / 36525, N = 125.04452 - 1934.136261jc
    (1.0004 - 0.0373cosd(N) + 0.0002cosd(2N)) * 0.69 * 0.2423 end
const M2_P = let T = (Dates.datetime2julian(start_date) - 2451545.0) / 36525, jc = (Dates.datetime2julian(DateTime(2015, 2, 15)) - 2451545.0) / 36525, N = 125.04452 - 1934.136261jc
    deg2rad(2(280.4664567 + 36000.76983T) - 2(218.3164477 + 481267.88123421T) - 2.14sind(N)) end
# M2_AMPLITUDE_SCALE (default 1) multiplies the M2_ONLY potential (diagnostics of rectified mean flows).
const M2_S = parse(Float64, get(ENV, "M2_AMPLITUDE_SCALE", "1"))
const M2_STOP = parse(Float64, get(ENV, "M2_FORCING_STOP_DAY", "Inf")) * 86400 # diagnostics: forcing off after this day
m2_potential(λ, φ, t) = -(t < M2_STOP) * M2_S * M2_G * M2_A * cosd(φ)^2 * cos(M2_ω * t + M2_P + 2 * deg2rad(λ))
Φ = m2_only ? m2_potential : Tidejinks.TidalPotential(grid, (start_date, end_date); reference_frame = "IAU_EARTH")

τ = parse(Float64, get(ENV, "BIHARMONIC_TIMESCALE", "0")) * days
@inline νhb(i, j, k, grid, τ) = Oceananigans.Operators.Azᶜᶜᶜ(i, j, k, grid)^2 / τ
ν = Field{Center, Center, Center}(grid)
τ > 0 && set!(ν, KernelFunctionOperation{Center, Center, Center}(νhb, grid, τ))
# BIHARMONIC_TIMESCALE (days; default 0 = off) and LAPLACIAN_VISCOSITY (m² s⁻¹, default 1000). With true depths, the biharmonic
# operator at steep partial-cell steps drives steady grid-scale η of tens to hundreds of metres, growing with ν.
νₗ = parse(Float64, get(ENV, "LAPLACIAN_VISCOSITY", "1000"))
closure = Tuple(c for c in (τ > 0 && get(ENV, "VISCOSITY", "default") != "biharmonic" ? HorizontalScalarBiharmonicDiffusivity(; ν) : nothing,
                            νₗ > 0 ? HorizontalScalarDiffusivity(ν = νₗ) : nothing) if !isnothing(c))
closure = isempty(closure) ? nothing : length(closure) == 1 ? only(closure) : closure
const Adapt = NumericalEarth.Oceans.Adapt
# VISCOSITY = mom6 replaces the above with Chen et al.'s MOM6 hor_visc settings: Laplacian KH = KH_VEL_SCALE Δ (0.03 m s⁻¹;
# RESOLN_SCALED_KH has no effect with one layer, where the baroclinic deformation radius is 0) and biharmonic
# AH = max(AH_VEL_SCALE Δ³ (0.01 m s⁻¹), SMAG_BI_CONST Δ⁴ |S| (0.06)), |S|² = (∂ₓu - ∂ᵧv)² + (∂ᵧu + ∂ₓv)², with
# Δ² = 2Δx²Δy² / (Δx² + Δy²) as in MOM6, both capped as by MOM6's BOUND_KH / BOUND_AH (KH ≤ Δ²/(8Δt), AH ≤ Δ⁴/(64Δt)).
# The Smagorinsky AH is recomputed at cell centres after every time step (lagged by one step). MOM6_VISCOSITY_SCALE
# (default 1) multiplies all three. VISCOSITY = flat uses the same operator with a uniform LAPLACIAN_VISCOSITY only.
# These are applied as a forcing with plain 2D (not thickness-weighted) componentwise operators, free slip at land
# (fluxes through land cells and land corners are zero, ∇²u is zeroed on boundary faces) and no metric terms:
# Oceananigans' closures weight both Laplacians of the biharmonic operator by the partial-cell thickness, so with true
# depths a 5 m cell next to a 4000 m one amplifies the inner Laplacian ~1000 times and the run blows up within ~10 steps
# (MOM6 forms the inner Laplacian without thickness weighting).
viscosity_name = get(ENV, "VISCOSITY", "default")
smagorinsky_update! = nothing
viscous_forcing = nothing
struct FlatViscosity{D, K, A}
    direction :: D
    kh :: K
    ah :: A
end
Adapt.@adapt_structure FlatViscosity
using Oceananigans.Grids: inactive_node, peripheral_node
using Oceananigans.Operators: Azᶠᶜᶜ, Azᶜᶠᶜ, Δxᶠᶜᶜ, Δyᶠᶜᶜ, Δxᶜᶠᶜ, Δyᶜᶠᶜ, ℑxyᶠᶠᵃ
const Cᵥ, Fᵥ = Center(), Face()
@inline wetc(i, j, grid) = !inactive_node(i, j, 1, grid, Cᵥ, Cᵥ, Cᵥ)
@inline wetu(i, j, grid) = !peripheral_node(i, j, 1, grid, Fᵥ, Cᵥ, Cᵥ)
@inline wetv(i, j, grid) = !peripheral_node(i, j, 1, grid, Cᵥ, Fᵥ, Cᵥ)
# Laplacian fluxes of a u-like field U (Face, Center) with coefficient ν (Center, Center field or nothing = 1)
@inline νc(i, j, ν) = @inbounds ν[i, j, 1]
@inline νc(i, j, ::Nothing) = 1
@inline νf(i, j, grid, ν) = ℑxyᶠᶠᵃ(i, j, 1, grid, ν)
@inline νf(i, j, grid, ::Nothing) = 1
@inline Fxu(i, j, grid, U, ν) = @inbounds ifelse(wetc(i, j, grid), νc(i, j, ν) * Δyᶜᶜᶜ(i, j, 1, grid) / Δxᶜᶜᶜ(i, j, 1, grid) * (U[i+1, j, 1] - U[i, j, 1]), 0.0)
@inline Fyu(i, j, grid, U, ν) = @inbounds ifelse(wetu(i, j, grid) & wetu(i, j-1, grid), νf(i, j, grid, ν) * Δxᶠᶠᶜ(i, j, 1, grid) / Δyᶠᶠᶜ(i, j, 1, grid) * (U[i, j, 1] - U[i, j-1, 1]), 0.0)
@inline Fxv(i, j, grid, V, ν) = @inbounds ifelse(wetv(i, j, grid) & wetv(i-1, j, grid), νf(i, j, grid, ν) * Δyᶠᶠᶜ(i, j, 1, grid) / Δxᶠᶠᶜ(i, j, 1, grid) * (V[i, j, 1] - V[i-1, j, 1]), 0.0)
@inline Fyv(i, j, grid, V, ν) = @inbounds ifelse(wetc(i, j, grid), νc(i, j, ν) * Δxᶜᶜᶜ(i, j, 1, grid) / Δyᶜᶜᶜ(i, j, 1, grid) * (V[i, j+1, 1] - V[i, j, 1]), 0.0)
@inline lapu(i, j, grid, U, ν) = ifelse(wetu(i, j, grid), (Fxu(i, j, grid, U, ν) - Fxu(i-1, j, grid, U, ν) + Fyu(i, j+1, grid, U, ν) - Fyu(i, j, grid, U, ν)) / Azᶠᶜᶜ(i, j, 1, grid), 0.0)
@inline lapv(i, j, grid, V, ν) = ifelse(wetv(i, j, grid), (Fxv(i+1, j, grid, V, ν) - Fxv(i, j, grid, V, ν) + Fyv(i, j, grid, V, ν) - Fyv(i, j-1, grid, V, ν)) / Azᶜᶠᶜ(i, j, 1, grid), 0.0)
# biharmonic: -∇·(A ∇L), L = ∇²U, fluxes of L as above
@inline BFxu(i, j, grid, U, A) = ifelse(wetc(i, j, grid), νc(i, j, A) * Δyᶜᶜᶜ(i, j, 1, grid) / Δxᶜᶜᶜ(i, j, 1, grid) * (lapu(i+1, j, grid, U, nothing) - lapu(i, j, grid, U, nothing)), 0.0)
@inline BFyu(i, j, grid, U, A) = ifelse(wetu(i, j, grid) & wetu(i, j-1, grid), νf(i, j, grid, A) * Δxᶠᶠᶜ(i, j, 1, grid) / Δyᶠᶠᶜ(i, j, 1, grid) * (lapu(i, j, grid, U, nothing) - lapu(i, j-1, grid, U, nothing)), 0.0)
@inline BFxv(i, j, grid, V, A) = ifelse(wetv(i, j, grid) & wetv(i-1, j, grid), νf(i, j, grid, A) * Δyᶠᶠᶜ(i, j, 1, grid) / Δxᶠᶠᶜ(i, j, 1, grid) * (lapv(i, j, grid, V, nothing) - lapv(i-1, j, grid, V, nothing)), 0.0)
@inline BFyv(i, j, grid, V, A) = ifelse(wetc(i, j, grid), νc(i, j, A) * Δxᶜᶜᶜ(i, j, 1, grid) / Δyᶜᶜᶜ(i, j, 1, grid) * (lapv(i, j+1, grid, V, nothing) - lapv(i, j, grid, V, nothing)), 0.0)
@inline bihu(i, j, grid, U, A) = ifelse(wetu(i, j, grid), (BFxu(i, j, grid, U, A) - BFxu(i-1, j, grid, U, A) + BFyu(i, j+1, grid, U, A) - BFyu(i, j, grid, U, A)) / Azᶠᶜᶜ(i, j, 1, grid), 0.0)
@inline bihv(i, j, grid, V, A) = ifelse(wetv(i, j, grid), (BFxv(i+1, j, grid, V, A) - BFxv(i, j, grid, V, A) + BFyv(i, j, grid, V, A) - BFyv(i, j-1, grid, V, A)) / Azᶜᶠᶜ(i, j, 1, grid), 0.0)
@inline bih_or_zero(f, i, j, grid, U, ::Nothing) = 0.0
@inline bih_or_zero(f, i, j, grid, U, A) = f(i, j, grid, U, A)
@inline lap_or_zero(f, i, j, grid, U, ::Nothing) = 0.0
@inline lap_or_zero(f, i, j, grid, U, ν) = f(i, j, grid, U, ν)
@inline (d::FlatViscosity{<:NumericalEarth.Oceans.XDirection})(i, j, k, grid, clock, fields) = lap_or_zero(lapu, i, j, grid, fields.u, d.kh) - bih_or_zero(bihu, i, j, grid, fields.u, d.ah)
@inline (d::FlatViscosity{<:NumericalEarth.Oceans.YDirection})(i, j, k, grid, clock, fields) = lap_or_zero(lapv, i, j, grid, fields.v, d.kh) - bih_or_zero(bihv, i, j, grid, fields.v, d.ah)
if viscosity_name ∈ ("mom6", "flat", "biharmonic")
    vscale = parse(Float64, get(ENV, "MOM6_VISCOSITY_SCALE", "1"))
    @inline Δ²(i, j, k, grid) = (dx = Δxᶜᶜᶜ(i, j, k, grid); dy = Δyᶜᶜᶜ(i, j, k, grid); 2dx^2 * dy^2 / (dx^2 + dy^2))
    @inline kh_mom6(i, j, k, grid, Δt, s) = min(s * 0.03 * sqrt(Δ²(i, j, k, grid)), Δ²(i, j, k, grid) / 8Δt)
    @inline kh_flat(i, j, k, grid, Δt, ν) = min(ν, Δ²(i, j, k, grid) / 8Δt)
    @inline ah_static(i, j, k, grid, Δt, s) = min(s * 0.01 * Δ²(i, j, k, grid)^1.5, Δ²(i, j, k, grid)^2 / 64Δt)
    # tension ∂ₓu - ∂ᵧv at centres; shear ∂ᵧu + ∂ₓv at corners, squared and averaged to centres
    @inline tension(i, j, k, grid, u, v) = δxᶜᶜᶜ(i, j, k, grid, u) / Δxᶜᶜᶜ(i, j, k, grid) - δyᶜᶜᶜ(i, j, k, grid, v) / Δyᶜᶜᶜ(i, j, k, grid)
    @inline shearᶠᶠᶜ(i, j, k, grid, u, v) = δyᶠᶠᶜ(i, j, k, grid, u) / Δyᶠᶠᶜ(i, j, k, grid) + δxᶠᶠᶜ(i, j, k, grid, v) / Δxᶠᶠᶜ(i, j, k, grid)
    @inline shear²ᶠᶠᶜ(i, j, k, grid, u, v) = shearᶠᶠᶜ(i, j, k, grid, u, v)^2
    @inline function ah_mom6(i, j, k, grid, Δt, s, u, v)
        S = sqrt(tension(i, j, k, grid, u, v)^2 + ℑxyᶜᶜᵃ(i, j, k, grid, shear²ᶠᶠᶜ, u, v))
        return min(s * max(0.01 * Δ²(i, j, k, grid)^1.5, 0.06 * Δ²(i, j, k, grid)^2 * S), Δ²(i, j, k, grid)^2 / 64Δt)
    end
    νₖ = viscosity_name == "biharmonic" ? nothing : Field{Center, Center, Center}(grid)
    viscosity_name == "biharmonic" || set!(νₖ, viscosity_name == "mom6" ? KernelFunctionOperation{Center, Center, Center}(kh_mom6, grid, Δt, vscale) :
                                        KernelFunctionOperation{Center, Center, Center}(kh_flat, grid, Δt, νₗ))
    isnothing(νₖ) || fill_halo_regions!(νₖ)
    νₐ = nothing
    if viscosity_name == "biharmonic"
        # VISCOSITY = biharmonic: A = Δ⁴/τ only (τ = BIHARMONIC_TIMESCALE days), bounded by Δ⁴/(64Δt)
        τ > 0 || error("VISCOSITY=biharmonic needs BIHARMONIC_TIMESCALE > 0 (days)")
        @inline ah_tau(i, j, k, grid, Δt, τ) = min(Δ²(i, j, k, grid)^2 / τ, Δ²(i, j, k, grid)^2 / 64Δt)
        νₐ = Field{Center, Center, Center}(grid); set!(νₐ, KernelFunctionOperation{Center, Center, Center}(ah_tau, grid, Δt, τ))
        fill_halo_regions!(νₐ)
        # Δ⁴/τ > Δ⁴/(64Δt) ⇔ τ < 64Δt, independent of Δ: the bound is hit everywhere or nowhere
        τ < 64Δt && @warn @sprintf("BIHARMONIC_TIMESCALE %g d < 64Δt = %g d: A bounded to Δ⁴/(64Δt) everywhere", τ / days, 64Δt / days)
        @info @sprintf("Flat biharmonic A = Δ⁴/τ, τ = %g days: A %.2e–%.2e m⁴ s⁻¹ (bounded: %s)", τ / days, extrema(interior(νₐ))..., τ < 64Δt)
    elseif viscosity_name == "mom6"
        νₐ = Field{Center, Center, Center}(grid); set!(νₐ, KernelFunctionOperation{Center, Center, Center}(ah_static, grid, Δt, vscale))
        fill_halo_regions!(νₐ)
        smagorinsky_update! = function (sim)
            u, v = sim.model.velocities.u, sim.model.velocities.v
            set!(νₐ, KernelFunctionOperation{Center, Center, Center}(ah_mom6, grid, Δt, vscale, u, v))
            fill_halo_regions!(νₐ)
            return nothing
        end
        @info @sprintf("MOM6 viscosity ×%g (flat operator): KH %.0f–%.0f m² s⁻¹, static AH %.2e–%.2e m⁴ s⁻¹", vscale,
                       extrema(interior(νₖ))..., extrema(interior(νₐ))...)
    elseif viscosity_name == "flat"
        @info @sprintf("Flat Laplacian viscosity: KH %.0f–%.0f m² s⁻¹", extrema(interior(νₖ))...)
    end
    closure = nothing
    viscous_forcing = (u = FlatViscosity(NumericalEarth.Oceans.XDirection(), νₖ, νₐ),
                       v = FlatViscosity(NumericalEarth.Oceans.YDirection(), νₖ, νₐ))
end
# FREE_SURFACE_KERNEL = shape (default, Shchepetkin & McWilliams), cosine or constant: the split-explicit averaging
# kernel; FREE_SURFACE_CFL (default 0.7) sets the substep count. Defaults reproduce NumericalEarth's default_free_surface.
fs_kernel = get(ENV, "FREE_SURFACE_KERNEL", "shape")
fs_cfl = parse(Float64, get(ENV, "FREE_SURFACE_CFL", "0.7"))
free_surface = if fs_kernel == "shape" && fs_cfl == 0.7
    NumericalEarth.Oceans.default_free_surface(grid; fixed_Δt = Δt)
else
    SE = Oceananigans.Models.HydrostaticFreeSurfaceModels.SplitExplicitFreeSurfaces
    kernel = fs_kernel == "cosine" ? SE.cosine_averaging_kernel : fs_kernel == "constant" ? SE.constant_averaging_kernel : SE.averaging_shape_function
    SplitExplicitFreeSurface(grid; cfl = fs_cfl, fixed_Δt = Δt, averaging_kernel = kernel)
end
@info "Free surface: $(summary(free_surface)), kernel $fs_kernel, cfl $fs_cfl"

loading = get(ENV, "SELF_ATTRACTION_AND_LOADING", "none")
loading_suffix = loading
self_attraction_and_loading = if loading == "scalar"
    ScalarSelfAttractionAndLoading(parse(Float64, get(ENV, "SCALAR_LOADING_COEFFICIENT", "0.1")))
elseif loading == "spherical_harmonic"
    @eval using SpeedyTransforms, ConservativeRegridding
    truncation = parse(Int, get(ENV, "LOADING_TRUNCATION", "40"))
    # LOADING_NLAT_HALF: Gaussian latitudes per hemisphere. LOADING_LOVE_NUMBERS = gegout (default), prem, or paper: the
    # kernel of CONVOLUTION_KERNEL=paper, (1 + k′ - h′ of the asymptotic fit) / 3, truncated like any other.
    nlat_half = parse(Int, get(ENV, "LOADING_NLAT_HALF", string(max(24, cld(truncation + 1, 2) + 1))))
    love_name = get(ENV, "LOADING_LOVE_NUMBERS", "gegout")
    love_numbers = if love_name == "paper"
        fit = NumericalEarth.Oceans.AsymptoticLoadLoveNumbers()
        (h = [1 - NumericalEarth.Oceans.one_plus_k_minus_h(fit, n) / 3 for n in 0:truncation], k = zeros(truncation + 1))
    elseif love_name == "paper_code"
        # The kernel as coded in Chen et al.'s MOM6 (signs differ from their Eq. 16), including its division by 3
        (h = [n == 0 ? 1.0 : 1 - (5.21196 + 8.865 / n) / 3 for n in 0:truncation], k = zeros(truncation + 1))
    elseif love_name == "prem"
        NumericalEarth.Oceans.prem_load_love_numbers()
    else
        nothing
    end
    SphericalHarmonicSelfAttractionAndLoading(grid; truncation, nlat_half, love_numbers)
elseif loading == "convolution"
    FT = get(ENV, "CONVOLUTION_PRECISION", "Float32") == "Float64" ? Float64 : Float32
    # CONVOLUTION_ACCURATE = true: remove_duplicates (the tripolar fold row's coincident cells) and subtract_uniform
    # (Float32 sums of differences, uniform-load response in Float64); see ConvolutionLoadingOperator.
    accurate = get(ENV, "CONVOLUTION_ACCURATE", "false") == "true"
    extra = accurate ? (; remove_duplicates = true, subtract_uniform = true) : (;)
    if get(ENV, "CONVOLUTION_KERNEL", "exact") == "paper"
        loading_suffix = "convolution_paper"
        ConvolutionSelfAttractionAndLoading(grid; FT, love_numbers = AsymptoticLoadLoveNumbers(), self_term = false,
                                            earth_density = 3 * 5517, extra...)
    else
        loading_suffix = "convolution_prem"
        ConvolutionSelfAttractionAndLoading(grid; FT, extra...)
    end
else
    nothing
end
@info "Self-attraction and loading: $(isnothing(self_attraction_and_loading) ? "none" : summary(self_attraction_and_loading))"

# Linear topographic wave drag -χ C / H u (Jayne & St. Laurent 2001) with WAVE_DRAG_COEFFICIENT = χ > 0, where
# C = κ h² N_b √(1 - f²/ω²) / 2 on 0.5ᵒ tiles comes from egbert_ray/solve/scripts/prep_wavedrag.py (zero where the
# tile is shallower than 100 m), interpolated to cell centres. The rate χ C / H is capped at 1 / 2Δt for stability.
struct WaveDrag{D, R}
    direction :: D
    rate :: R
end
Adapt.@adapt_structure WaveDrag
@inline (d::WaveDrag{<:NumericalEarth.Oceans.XDirection})(i, j, k, grid, clock, fields) = @inbounds - (d.rate[i-1, j, 1] + d.rate[i, j, 1]) / 2 * fields.u[i, j, k]
@inline (d::WaveDrag{<:NumericalEarth.Oceans.YDirection})(i, j, k, grid, clock, fields) = @inbounds - (d.rate[i, j-1, 1] + d.rate[i, j, 1]) / 2 * fields.v[i, j, k]

#
# WAVE_DRAG_FIELD = jsl_rh: C is MOM6's JSL_RH (Chen et al., JSLH_p04arc9_aw.nc, regridded; needs EXTERNAL_DEPTH = mom6 for
# the file) and χ plays BT_WAVE_DRAG_SCALE. WAVE_DRAG_IMPLICIT = true: no cap, and the rate r = χ C / H is replaced by
# r / (1 + r Δt), so that a forward step gives u (1 - r Δt / (1 + r Δt)) = u / (1 + r Δt), the backward-Euler (semi-implicit)
# damping, which is stable for any r (in the AB2/RK sense the effective rate never exceeds 1/Δt).
χ = parse(Float64, get(ENV, "WAVE_DRAG_COEFFICIENT", "0"))
wave_drag_field = get(ENV, "WAVE_DRAG_FIELD", "egbert_ray")
wave_drag_implicit = get(ENV, "WAVE_DRAG_IMPLICIT", "false") == "true"
forcing = if χ > 0
    C = reshape(reinterpret(Float64, read("wavedrag_halfdeg_C.bin")), 720, 312) # 0.5ᵒ tiles, lat from -77.75 to 77.75
    function tile_value(λ, φ)
        i = clamp(floor(Int, mod(λ, 360) / 0.5) + 1, 1, 720)
        j = floor(Int, (φ + 78) / 0.5) + 1
        return 1 ≤ j ≤ 312 ? C[i, j] : 0.0
    end
    rh = wave_drag_field == "jsl_rh" ? load(@sprintf("mom6_p04arc9_on_%d.jld2", Nx), "rh") : nothing
    rate = zeros(Nx, Ny)
    raw_rate = zeros(Nx, Ny)
    for j in 1:Ny, i in 1:Nx
        H = -bottom_height[i, j, 1]
        λ = Oceananigans.Grids.λnode(i, j, 1, cpu_grid, Center(), Center(), Center())
        φ = Oceananigans.Grids.φnode(i, j, 1, cpu_grid, Center(), Center(), Center())
        c = isnothing(rh) ? tile_value(λ, φ) : rh[i, j]
        raw_rate[i, j] = H > 0 ? χ * c / H : 0.0
        rate[i, j] = wave_drag_implicit ? raw_rate[i, j] / (1 + raw_rate[i, j] * Δt) : min(raw_rate[i, j], 1 / 2Δt)
    end
    @info @sprintf("Wave drag %s χ = %g: max r = %.3e s⁻¹ (1/r = %.2f h), max r Δt = %.3g, %d cells with r Δt > 0.5, implicit = %s",
                   wave_drag_field, χ, maximum(raw_rate), 1 / maximum(raw_rate) / 3600, maximum(raw_rate) * Δt,
                   count(raw_rate .* Δt .> 0.5), wave_drag_implicit)
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
if !isnothing(viscous_forcing)
    forcing = isempty(forcing) ? viscous_forcing : (u = (forcing.u, viscous_forcing.u), v = (forcing.v, viscous_forcing.v))
end
wave_drag_suffix = χ > 0 ? @sprintf("_wavedrag%s%g", wave_drag_field == "jsl_rh" ? "jsl" : "", χ) : ""

# CORIOLIS_SCHEME = EnstrophyConserving (default), EnergyConserving, ActiveWeightedEnstrophyConserving,
# ActiveWeightedEnergyConserving or TriadScheme
# or ThicknessWeightedEnergy: Oceananigans' EnergyConserving form -ℑʸ(f ℑˣ(A_y v)) / A_y^u conserves energy only when
# the cell thickness Δz is uniform; with partial cells of very different thickness (5 m next to 4000 m) it does work,
# and a steady geostrophic mode grows near steep Antarctic slopes without any forcing. ThicknessWeightedEnergy uses
# the Sadourny energy form on transports with potential vorticity f / h at corners,
#   Fᵘ = -ℑʸ(f ℑˣ(A_y v) / hᶠᶠ) A_x^u / V^u,   Fᵛ = +ℑˣ(f ℑʸ(A_x u) / hᶠᶠ) A_y^v / V^v,
# hᶠᶠ the mean thickness of the active cells around the corner, whose work cancels exactly, applied as a forcing.
scheme_name = get(ENV, "CORIOLIS_SCHEME", "EnstrophyConserving")
thickness_weighted_coriolis = scheme_name == "ThicknessWeightedEnergy"
coriolis_scheme = thickness_weighted_coriolis ? nothing :
    getproperty(scheme_name ∈ ("EnstrophyConserving", "EnergyConserving") ? Oceananigans.Advection : Oceananigans.Coriolis, Symbol(scheme_name))()
coriolis = get(ENV, "CORIOLIS", "on") == "off" || thickness_weighted_coriolis ? nothing : HydrostaticSphericalCoriolis(scheme = coriolis_scheme) # CORIOLIS=off: diagnostics only
if thickness_weighted_coriolis
    struct ThicknessWeightedCoriolis{D, F}
        direction :: D
        f :: F # Coriolis parameter at (Face, Face)
    end
    Adapt.@adapt_structure ThicknessWeightedCoriolis
    @inline Ayvᶜᶠᶜ(i, j, k, grid, v) = @inbounds ifelse(wetv(i, j, grid), Ayᶜᶠᶜ(i, j, k, grid) * v[i, j, k], 0.0)
    @inline Axuᶠᶜᶜ(i, j, k, grid, u) = @inbounds ifelse(wetu(i, j, grid), Axᶠᶜᶜ(i, j, k, grid) * u[i, j, k], 0.0)
    @inline wetΔz(i, j, k, grid) = ifelse(wetc(i, j, grid), Δzᶜᶜᶜ(i, j, k, grid), 0.0)
    @inline wet1(i, j, k, grid) = ifelse(wetc(i, j, grid), 1.0, 0.0)
    # CORIOLIS_CORNER_THICKNESS = faces (default): h is the mean thickness of the wet velocity faces meeting at the corner
    # (consistent with the face thicknesses in A_x, A_y, so f ℑˣ(A_y v)/h ≈ f ℑˣ(Δx v) on slopes); cells: mean of the
    # wet cells around the corner (biased low against the min-of-neighbours face thickness on steep slopes).
    @inline wetΔzu(i, j, k, grid) = ifelse(wetu(i, j, grid), Oceananigans.Operators.Δzᶠᶜᶜ(i, j, k, grid), 0.0)
    @inline wetΔzv(i, j, k, grid) = ifelse(wetv(i, j, grid), Oceananigans.Operators.Δzᶜᶠᶜ(i, j, k, grid), 0.0)
    @inline wet1u(i, j, grid) = ifelse(wetu(i, j, grid), 1.0, 0.0)
    @inline wet1v(i, j, grid) = ifelse(wetv(i, j, grid), 1.0, 0.0)
    @inline function hᶠᶠ_faces(i, j, k, grid)
        n = wet1u(i, j, grid) + wet1u(i, j-1, grid) + wet1v(i, j, grid) + wet1v(i-1, j, grid)
        s = wetΔzu(i, j, k, grid) + wetΔzu(i, j-1, k, grid) + wetΔzv(i, j, k, grid) + wetΔzv(i-1, j, k, grid)
        return ifelse(n > 0, s / max(n, 1), Inf)
    end
    @inline function hᶠᶠ_cells(i, j, k, grid)
        n = wet1(i, j, k, grid) + wet1(i-1, j, k, grid) + wet1(i, j-1, k, grid) + wet1(i-1, j-1, k, grid)
        s = wetΔz(i, j, k, grid) + wetΔz(i-1, j, k, grid) + wetΔz(i, j-1, k, grid) + wetΔz(i-1, j-1, k, grid)
        return ifelse(n > 0, s / max(n, 1), Inf)
    end
    if get(ENV, "CORIOLIS_CORNER_THICKNESS", "faces") == "cells"
        @inline hᶠᶠ(i, j, k, grid) = hᶠᶠ_cells(i, j, k, grid)
    else
        @inline hᶠᶠ(i, j, k, grid) = hᶠᶠ_faces(i, j, k, grid)
    end
    @inline qVᶠᶠ(i, j, k, grid, f, v) = @inbounds f[i, j, 1] * ℑxᶠᵃᵃ(i, j, k, grid, Ayvᶜᶠᶜ, v) / hᶠᶠ(i, j, k, grid)
    @inline qUᶠᶠ(i, j, k, grid, f, u) = @inbounds f[i, j, 1] * ℑyᵃᶠᵃ(i, j, k, grid, Axuᶠᶜᶜ, u) / hᶠᶠ(i, j, k, grid)
    @inline (c::ThicknessWeightedCoriolis{<:NumericalEarth.Oceans.XDirection})(i, j, k, grid, clock, fields) =
        ifelse(wetu(i, j, grid), ℑyᵃᶜᵃ(i, j, k, grid, qVᶠᶠ, c.f, fields.v) * Axᶠᶜᶜ(i, j, k, grid) / Vᶠᶜᶜ(i, j, k, grid), 0.0)
    @inline (c::ThicknessWeightedCoriolis{<:NumericalEarth.Oceans.YDirection})(i, j, k, grid, clock, fields) =
        ifelse(wetv(i, j, grid), - ℑxᶜᵃᵃ(i, j, k, grid, qUᶠᶠ, c.f, fields.u) * Ayᶜᶠᶜ(i, j, k, grid) / Vᶜᶠᶜ(i, j, k, grid), 0.0)
    fᶠᶠ = Field{Face, Face, Nothing}(grid)
    set!(fᶠᶠ, (λ, φ) -> 2 * 7.292115e-5 * sind(φ))
    fill_halo_regions!(fᶠᶠ)
    coriolis_forcing = (u = ThicknessWeightedCoriolis(NumericalEarth.Oceans.XDirection(), fᶠᶠ),
                        v = ThicknessWeightedCoriolis(NumericalEarth.Oceans.YDirection(), fᶠᶠ))
    astuple(x) = x isa Tuple ? x : (x,)
    forcing = isempty(forcing) ? coriolis_forcing :
              (u = (astuple(forcing.u)..., coriolis_forcing.u), v = (astuple(forcing.v)..., coriolis_forcing.v))
end
coriolis_suffix = scheme_name == "EnstrophyConserving" ? "" : "_" * scheme_name

bottom_drag_coefficient = parse(Float64, get(ENV, "BOTTOM_DRAG_COEFFICIENT", "0.003"))
# MOMENTUM_ADVECTION = none (default: linear momentum, as WENO vector-invariant advection blows up next to 20 m cells) or default
momentum_advection = get(ENV, "MOMENTUM_ADVECTION", "none") == "none" ? (; momentum_advection = nothing) : (;)
# BOTTOM_DRAG_FORM = imex (default: NumericalEarth's IMEX flux BC, whose drag is carried entirely by the vertically
# implicit solve) or explicit (FluxBoundaryCondition, enters the slow tendency and hence the barotropic mode).
drag_form = get(ENV, "BOTTOM_DRAG_FORM", "imex")
ocean = ocean_simulation(grid; momentum_advection..., Δt, implicit_bottom_drag = drag_form == "imex", free_surface, barotropic_potential = Φ, closure, bottom_drag_coefficient, forcing, coriolis,
                         self_attraction_and_loading, equation_of_state = LinearEquationOfState(),
                         radiative_forcing = nothing, gravitational_acceleration = g, stop_time = stop_days * days)
set!(ocean.model, T = 10, S = 35)
isnothing(smagorinsky_update!) || add_callback!(ocean, smagorinsky_update!, IterationInterval(1))
# SEED_VELOCITY = a (m/s, diagnostics): start from u, v = a × white noise on wet faces (seeded RNG).
# ENERGY_DIAGNOSTIC = true: log E = ½ Σ (Vᵘ u² + Vᵛ v²) + ½ g Σ A_z η² (per unit density) at every progress interval.
seed = parse(Float64, get(ENV, "SEED_VELOCITY", "0"))
if seed > 0
    let rng = Base.require(Base.PkgId(Base.UUID("9a3f8284-a2c9-5f02-9a11-845980a1fd5c"), "Random"))
        rng.seed!(1)
        u, v = ocean.model.velocities.u, ocean.model.velocities.v
        set!(u, seed .* randn(size(u)...)); set!(v, seed .* randn(size(v)...))
        Oceananigans.ImmersedBoundaries.mask_immersed_field!(u); Oceananigans.ImmersedBoundaries.mask_immersed_field!(v)
        fill_halo_regions!((u, v))
    end
end
energy_diagnostic = get(ENV, "ENERGY_DIAGNOSTIC", "false") == "true"
if energy_diagnostic
    Vu = Array(interior(Field(KernelFunctionOperation{Face, Center, Center}(Vᶠᶜᶜ, grid))))
    Vv = Array(interior(Field(KernelFunctionOperation{Center, Face, Center}(Vᶜᶠᶜ, grid))))
    Az = Array(interior(Field(KernelFunctionOperation{Center, Center, Center}(Oceananigans.Operators.Azᶜᶜᶜ, grid))))[:, :, 1]
    wet_c = Array(interior(Field(KernelFunctionOperation{Center, Center, Center}((i, j, k, g) -> wetc(i, j, g) ? 1.0 : 0.0, grid))))[:, :, 1]
    function total_energy(sim)
        u = Array(interior(sim.model.velocities.u)); v = Array(interior(sim.model.velocities.v)); e = Array(interior(η, :, :, 1))
        K = sum(ifelse.(isfinite.(u), Vu .* u .^ 2, 0.0)) / 2 + sum(ifelse.(isfinite.(v), Vv .* v .^ 2, 0.0)) / 2
        P = g * sum(Az .* wet_c .* e .^ 2) / 2
        return K, P
    end
    add_callback!(ocean, sim -> (KP = total_energy(sim); @info @sprintf("ENERGY day %.3f: K %.6e P %.6e E %.6e m⁵ s⁻²", time(sim) / 86400, KP..., sum(KP))),
                  TimeInterval(parse(Float64, get(ENV, "ENERGY_INTERVAL_DAYS", "0.25")) * days))
end
# SAL_UPDATE_INTERVAL = N > 1 (time steps, default 1): recompute the SAL potential or gradient only every N steps and hold it
# in between (for the O(N²) convolution, which costs ~7 s per update in Float32 on a T4 at 0.36ᵒ).
sal_interval = parse(Int, get(ENV, "SAL_UPDATE_INTERVAL", "1"))
if haskey(ENV, "SAL_UPDATE_INTERVAL") && haskey(ocean.callbacks, :self_attraction_and_loading) # (=1: once per step rather than per substage)
    # once per N steps at the end of a step (the UpdateStateCallsite callback runs at every substage of the time stepper)
    delete!(ocean.callbacks, :self_attraction_and_loading)
    add_callback!(ocean, sim -> NumericalEarth.Oceans.compute_self_attraction_and_loading!(sim.model, self_attraction_and_loading),
                  IterationInterval(sal_interval); name = :self_attraction_and_loading_held)
    @info "SAL updated every $sal_interval steps ($(sal_interval * Δt / 60) min)"
end

η = ocean.model.free_surface.displacement
# RUN_TAG is appended to the output name.
name = @sprintf("barotropic_tides_%gdeg_dt%gmin_%s%s%s", 1 / cells_per_degree, Δt / minutes, loading_suffix,
                wave_drag_suffix * coriolis_suffix, get(ENV, "RUN_TAG", ""))
# FIT_START_DAY = d > 0: instead of hourly η, accumulate a least-squares fit of a mean and M2, N2, S2, K1, O1 to hourly
# η from day d to the end (N2 separates from M2 only over ≥ 27.6 days), and save the complex amplitudes
# Z (η ≈ Re Z exp(iωt), t from start_date) in output/<name>_harmonics.jld2, with η every 5 days for checking.
fit_start = parse(Float64, get(ENV, "FIT_START_DAY", "0")) * days
fit_periods = m2_only ? (; M2 = 12.4206012) : (M2 = 12.4206012, N2 = 12.65834751, S2 = 12.0, K1 = 23.9344696, O1 = 25.8193417) # hours
fit_ω = [2π / (3600T) for T in values(fit_periods)]
fit_times = Float64[]
fit_sums = [Field{Center, Center, Nothing}(grid) for _ in 1:(2length(fit_ω) + 1)]
# FIT_VELOCITIES=true also fits u and v (saved as Zu, Zv, with Oceananigans' Face/Center interiors) for energy budgets.
fit_velocities = get(ENV, "FIT_VELOCITIES", "false") == "true"
fit_sums_u = [Field{Face, Center, Center}(grid) for _ in 1:(fit_velocities ? 2length(fit_ω) + 1 : 0)]
fit_sums_v = [Field{Center, Face, Center}(grid) for _ in 1:(fit_velocities ? 2length(fit_ω) + 1 : 0)]
η_snapshots = Dict{String, Array{Float32}}()
function accumulate_fit(sim)
    t = time(sim)
    t < fit_start - 1 && return nothing
    push!(fit_times, t)
    basis = [1.0; cos.(fit_ω .* t); sin.(fit_ω .* t)]
    for (S, b) in zip(fit_sums, basis)
        interior(S) .+= b .* interior(η)
    end
    if fit_velocities
        u, v = sim.model.velocities.u, sim.model.velocities.v
        for (Su, Sv, b) in zip(fit_sums_u, fit_sums_v, basis)
            interior(Su) .+= b .* interior(u)
            interior(Sv) .+= b .* interior(v)
        end
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
    extra = (;)
    if fit_velocities
        fitted(sums) = begin
            F = [Array(interior(f, :, :, 1)) for f in sums]
            c(n) = sum(G⁻¹[n, m] .* F[m] for m in eachindex(F))
            Dict(string(q) => ComplexF32.(c(1 + k) .- im .* c(1 + K + k)) for (k, q) in enumerate(keys(fit_periods)))
        end
        meanof(sums) = (F = [Array(interior(f, :, :, 1)) for f in sums]; Float32.(sum(G⁻¹[1, m] .* F[m] for m in eachindex(F))))
        extra = (; Zu = fitted(fit_sums_u), Zv = fitted(fit_sums_v), mean_u = meanof(fit_sums_u), mean_v = meanof(fit_sums_v))
    end
    jldsave(joinpath("output", name * "_harmonics.jld2"); Z, mean = coefficient(1), first_day = t[1] / 86400,
            last_day = t[end] / 86400, samples = length(t), η_snapshots, extra...)
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
    time(sim) > 1.5days && push!(day_walls, (elapsed, iteration(sim)))
    flush(stderr)
    wall_time[] = time_ns()
end
day_walls = Tuple{Float64, Int}[]
add_callback!(ocean, progress, TimeInterval(parse(Float64, get(ENV, "PROGRESS_INTERVAL_DAYS", "1")) * days))

run!(ocean)
fit_start > 0 && save_fit(ocean)

# Cost: wall seconds per model day and per step from the progress intervals after day 1 (excluding compilation), and the
# SAL update alone (compute_self_attraction_and_loading!, synchronised), averaged over 20 calls after one warm-up.
if length(day_walls) ≥ 2
    w = sum(first, day_walls[2:end]); n = day_walls[end][2] - day_walls[1][2]
    interval = parse(Float64, get(ENV, "PROGRESS_INTERVAL_DAYS", "1"))
    @info @sprintf("COST %s: %.3f wall s per model day, %.2f ms per step (%d steps)", name, w / (interval * (length(day_walls) - 1)), 1e3w / n, n)
end
sync() = arch_name == "CUDA" ? CUDA.synchronize() : nothing
if !isnothing(self_attraction_and_loading) && hasmethod(NumericalEarth.Oceans.compute_self_attraction_and_loading!,
                                                        Tuple{Any, typeof(self_attraction_and_loading)})
    f! = NumericalEarth.Oceans.compute_self_attraction_and_loading!
    f!(ocean.model, self_attraction_and_loading); sync()
    t0 = time_ns(); for _ in 1:20; f!(ocean.model, self_attraction_and_loading); end; sync()
    @info @sprintf("SALCOST %s: %.2f ms per SAL update", name, (time_ns() - t0) / 20 / 1e6)
end
