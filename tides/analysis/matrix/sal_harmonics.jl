# Offline SAL of the M2 harmonic: apply the run's SphericalHarmonicSelfAttractionAndLoading operator (truncation 40,
# nlat_half 24, PREM load Love numbers, on the run's immersed grid: MOM6 aw depth, 5 m minimum, 9314 m layer) to
# pᵇ = g Re Z and pᵇ = g Im Z. In the model pᵇ = g η − Σ b Δz minus its initial value; T, S are uniform and fixed, so
# pᵇ = g η exactly. Ψ is linear in pᵇ, so Ψ_M2 = Ψ[g Re Z] + i Ψ[g Im Z], same time convention as Z (η = Re Z e^{iωt}).
using NumericalEarth, Oceananigans, SpeedyTransforms, ConservativeRegridding, JLD2, Printf
using Oceananigans.BoundaryConditions: fill_halo_regions!
const O = NumericalEarth.Oceans
root = joinpath(@__DIR__, "..", "..", "GlobalKiloMischief", "tides")
name = "barotropic_tides_0.36deg_dt5min_spherical_harmonic_wavedragjsl0.8_ThicknessWeightedEnergy_cal_tau30_hf"
g = 9.8; Nx, Ny = 1000, 472; minimum_depth = 5.0
D = load(joinpath(root, "mom6_p04arc9_on_1000.jld2"), "depth")
Hmax = max(6000.0, ceil(maximum(D)))
ug = TripolarGrid(CPU(); size = (Nx, Ny, 1), halo = (7, 7, 4), z = (-Hmax, 0))
zb = Field{Center, Center, Nothing}(ug); set!(zb, reshape(map(x -> x > 0 ? -max(x, minimum_depth) : 0.0, D), Nx, Ny, 1)); fill_halo_regions!(zb)
grid = ImmersedBoundaryGrid(ug, PartialCellBottom(zb; minimum_fractional_cell_height = min(1e-3, minimum_depth / 2Hmax)); active_cells_map = true)
truncation = 40; nlat_half = max(24, cld(truncation + 1, 2) + 1)
sal = SphericalHarmonicSelfAttractionAndLoading(grid; truncation, nlat_half, love_numbers = O.prem_load_love_numbers())
Z = load(joinpath(root, "output", name * "_harmonics.jld2"), "Z")["M2"]
wet = D .> 0
function Ψof(η)
    p = sal.bottom_pressure_anomaly
    set!(p, reshape(g .* ifelse.(wet, η, 0.0), Nx, Ny, 1)); fill_halo_regions!(p)
    O.apply_loading_operator!(sal.potential, sal.operator, p)
    Array(interior(sal.potential, :, :, 1))
end
Psi = complex.(Ψof(real.(Z)), Ψof(imag.(Z))) ./ g
Psi[.!wet] .= 0
area, φ, Hours = load(joinpath(@__DIR__, "..", "diag", "grid_1000.jld2"), "area", "φ", "depth")
deep = wet .& (D .> 1000) .& (abs.(φ) .< 66)
r(x) = sqrt(sum(abs2.(x[deep]) .* area[deep]) / sum(area[deep]))
proj = sum((Psi .* conj.(Z) .* area)[deep]) / sum((abs2.(Z) .* area)[deep])
@printf("deep rms|Ψ/g| / rms|Z| = %.4f; projection Ψ/g on Z = %.4f %+.4fi (≈ −β for a scalar approximation)\n", r(Psi) / r(Z), real(proj), imag(proj))
out = joinpath(root, "output", name * "_sal_harmonics.jld2")
jldsave(out; Psi_M2 = ComplexF64.(Psi), Z_M2 = ComplexF64.(Z),
        settings = (operator = "SphericalHarmonicSelfAttractionAndLoading", truncation, nlat_half, love_numbers = "PREM (Wang et al. 2012) via prem_load_love_numbers()",
                    g, pb = "kinematic pᵇ = g η (uniform T, S; initial-state reference removed)",
                    units = "metres: Psi_M2 = Ψ/g",
                    sign = "Ψ is added to the tidal potential; momentum forcing is −∇Ψ = −g ∇(Psi_M2). Ψ = −Σₗ γₗ pᵇₗ, so Psi_M2 ≈ −β Z with β ≈ 0.1",
                    time = "η = Re Z e^{iωt}, Ψ = Re Ψ_M2 e^{iωt}, t from 2015-01-01 00:00",
                    grid = "TripolarGrid 1000×472, MOM6 aw depth min 5 m, layer $(Hmax) m, PartialCellBottom ϵ = $(min(1e-3, minimum_depth / 2Hmax))"))
println("saved ", out)
