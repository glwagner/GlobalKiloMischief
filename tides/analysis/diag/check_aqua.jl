# Analytic verification: non-rotating aquaplanet of uniform depth, forced by Φ = A P_n²(sinφ) cos(2λ - ωt) with uniform linear
# drag r. For degree n: Z = -HκA / (g-free) ... η̂ = -Hκ A / (-ω² + iωr + gHκ(1-γ_n)), κ = n(n+1)/R². Any amplitude or phase
# departure measures forcing, SAL, split-explicit and time-step errors, including a time lag in SAL.
# julia --project=sal_env check_aqua.jl <n> <SAL> <Δt min> [substeps]
using NumericalEarth, Oceananigans, Oceananigans.Units, SpeedyTransforms, ConservativeRegridding, CUDA, Printf
using Oceananigans.Grids: λnode, φnode
using Oceananigans.Operators: Azᶜᶜᶜ
const O = NumericalEarth.Oceans; const Adapt = O.Adapt
n = parse(Int, ARGS[1]); sal_name = ARGS[2]; Δt = parse(Float64, ARGS[3]) * 60
const ω = 2π / (12.4206012 * 3600); const A₀ = 1.0; R = 6371e3; Hd = 4000.0; r = 1 / (3 * 86400.0); g = 9.80665
# Angular factors: n = 2: cos²φ; n = 4: cos²φ (7 sin²φ - 1)
const deg4 = n == 4
angular(φ) = deg4 ? cosd(φ)^2 * (7sind(φ)^2 - 1) : cosd(φ)^2
Φf(λ, φ, t) = ifelse(isfinite(λ + φ), A₀ * angular(φ) * cos(2 * deg2rad(λ) - ω * t), zero(φ))
struct LinDrag{D, FT}; d::D; r::FT; end
Adapt.@adapt_structure LinDrag
@inline (f::LinDrag{<:O.XDirection})(i, j, k, grid, clock, fields) = @inbounds -f.r * fields.u[i, j, k]
@inline (f::LinDrag{<:O.YDirection})(i, j, k, grid, clock, fields) = @inbounds -f.r * fields.v[i, j, k]
Nx, Ny = 360, 170
ug = TripolarGrid(GPU(); size = (Nx, Ny, 1), halo = (7, 7, 4), z = (-6000, 0))
zb = Field{Center, Center, Nothing}(ug); cg0 = TripolarGrid(CPU(); size = (Nx, Ny, 1), halo = (7, 7, 4), z = (-6000, 0))
φ0 = [φnode(i, j, 1, cg0, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
bed = [(j ≤ 2 || j ≥ Ny - 1 || !isfinite(φ0[i, j]) || abs(φ0[i, j]) > 78) ? 0.0 : -Hd for i in 1:Nx, j in 1:Ny]
set!(zb, reshape(bed, Nx, Ny, 1)); Oceananigans.BoundaryConditions.fill_halo_regions!(zb)
grid = ImmersedBoundaryGrid(ug, PartialCellBottom(zb; minimum_fractional_cell_height = 1e-3); active_cells_map = true)
sal = sal_name == "none" ? nothing : sal_name == "scalar" ? ScalarSelfAttractionAndLoading(0.1) :
      sal_name == "sh" ? SphericalHarmonicSelfAttractionAndLoading(grid; truncation = 40, nlat_half = 24) :
      ConvolutionSelfAttractionAndLoading(grid; FT = Float32)
free_surface = length(ARGS) > 3 ? SplitExplicitFreeSurface(grid; substeps = parse(Int, ARGS[4])) : O.default_free_surface(grid; fixed_Δt = Δt)
days_run = 16.0
ocean = ocean_simulation(grid; momentum_advection = nothing, Δt, free_surface, barotropic_potential = Φf, closure = HorizontalScalarDiffusivity(ν = 1.0),
                         bottom_drag_coefficient = 1e-7, coriolis = nothing, self_attraction_and_loading = sal,
                         forcing = (u = LinDrag(O.XDirection(), r), v = LinDrag(O.YDirection(), r)), equation_of_state = LinearEquationOfState(),
                         radiative_forcing = nothing, gravitational_acceleration = g, stop_time = days_run * days, verbose = false)
set!(ocean.model, T = 10, S = 35)
η = ocean.model.free_surface.displacement
a = zeros(Nx, Ny); b = zeros(Nx, Ny); N = Ref(0)
function acc(sim)
    t = time(sim); t < 8days && return
    e = Array(interior(η, :, :, 1)); a .+= e .* cos(ω * t); b .+= e .* sin(ω * t); N[] += 1
end
add_callback!(ocean, acc, TimeInterval(1hours))
add_callback!(ocean, s -> (e = Array(interior(η, :, :, 1)); println("day ", time(s) / 86400, " maxη ", maximum(abs, filter(isfinite, e)), " nan ", count(isnan, e)); flush(stdout)), TimeInterval(1days))
run!(ocean)
Z = 2 .* (a .- im .* b) ./ N[]
cg = TripolarGrid(CPU(); size = (Nx, Ny, 1), halo = (7, 7, 4), z = (-6000, 0))
λ = [λnode(i, j, 1, cg, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]; φ = [φnode(i, j, 1, cg, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
ar = [Azᶜᶜᶜ(i, j, 1, cg) for i in 1:Nx, j in 1:Ny]
Y = @. angular(φ) * cis(2 * deg2rad(λ))
good = isfinite.(Z) .& isfinite.(ar) .& isfinite.(Y) .& (abs.(φ) .< 70); println("good ", count(good), " of ", length(Z), " finite Z ", count(isfinite, Z)); ar = ar .* good; Z = ifelse.(good, Z, 0); Y = ifelse.(good, Y, 0)
# Force sign: ∂t u = -∇(gη(1-γ) + Φ); Φ̂ e^{iωt} with Φ̂ = A Y*, i.e. Re Ŷ e^{-i(…)}: match by projecting on conj(Y)
love = Dict(2 => 1.686, 4 => 2 - 0.1337 + 0.0534 + 0.0) # placeholder replaced below
lines = [split(l) for l in eachline("/shared/home/greg/tides/egbert_ray/data/LoveNumbers/gegout-lln.txt") if !startswith(strip(l), '#')]
Lg = 1 + parse(Float64, lines[n+1][3]) - parse(Float64, lines[n+1][2]); Lw = Lg
γ = sal_name == "none" ? 0.0 : sal_name == "scalar" ? 0.1 : 3 * 1035 / 5517 * Lg / (2n + 1)
κ = n * (n + 1) / R^2
exact(γ) = -Hd * κ * A₀ / (-ω^2 + im * ω * r + g * Hd * κ * (1 - γ))
proj = sum(conj.(conj.(Y)) .* Z .* ar) / sum(abs2.(Y) .* ar)  # coefficient of conj(Y) e^{iωt}
res = sqrt(sum(abs2.(Z .- proj .* conj.(Y)) .* ar) / sum(abs2.(proj .* conj.(Y)) .* ar))
e = exact(γ)
@printf("RESULT n=%d SAL=%s Δt=%.1f min: |Z|=%.5f exact %.5f ratio %.4f; phase %.3f exact %.3f diff %.3f deg; residual %.4f; lag-equivalent %.1f s\n", n, sal_name, Δt / 60, abs(proj), abs(e), abs(proj) / abs(e), rad2deg(angle(proj)), rad2deg(angle(e)), rad2deg(angle(proj / e)), res, -angle(proj / e) / ω)
