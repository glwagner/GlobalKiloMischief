# Work done by the discrete Coriolis term on a random velocity field, native depth with 5 m minimum (no cell removal):
# W = Σ ρ (Vᵘ u Fᵘ + Vᵛ v Fᵛ) for Oceananigans' EnergyConserving and for ThicknessWeightedEnergy (as in barotropic_tides.jl).
# Saves per-cell work maps (u and v contributions summed onto cells) to coriolis_work.jld2.
using NumericalEarth, Oceananigans, JLD2, Printf, Random
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Operators
using Oceananigans.Grids: peripheral_node, inactive_node
using Oceananigans.Coriolis: x_f_cross_U, y_f_cross_U
using Oceananigans.Advection: EnergyConserving
Nx, Ny = 1000, 472
ug = TripolarGrid(CPU(); size = (Nx, Ny, 1), halo = (7, 7, 4), z = (-6000, 0))
bh = regrid_bathymetry(ug; minimum_depth = 5, interpolation_passes = 3, major_basins = 1)
grid = ImmersedBoundaryGrid(ug, PartialCellBottom(bh; minimum_fractional_cell_height = 5 / 12000); active_cells_map = true)
const C, F = Center(), Face()
wetc(i, j, g) = !inactive_node(i, j, 1, g, C, C, C)
wetu(i, j, g) = !peripheral_node(i, j, 1, g, F, C, C)
wetv(i, j, g) = !peripheral_node(i, j, 1, g, C, F, C)
fᶠᶠ = Field{Face, Face, Nothing}(grid); set!(fᶠᶠ, (λ, φ) -> 2 * 7.292115e-5 * sind(φ)); fill_halo_regions!(fᶠᶠ)
Ayv(i, j, k, g, v) = ifelse(wetv(i, j, g), Ayᶜᶠᶜ(i, j, k, g) * v[i, j, k], 0.0)
Axu(i, j, k, g, u) = ifelse(wetu(i, j, g), Axᶠᶜᶜ(i, j, k, g) * u[i, j, k], 0.0)
wΔz(i, j, k, g) = ifelse(wetc(i, j, g), Δzᶜᶜᶜ(i, j, k, g), 0.0); w1(i, j, k, g) = ifelse(wetc(i, j, g), 1.0, 0.0)
hff(i, j, k, g) = (n = w1(i, j, k, g) + w1(i-1, j, k, g) + w1(i, j-1, k, g) + w1(i-1, j-1, k, g);
                   s = wΔz(i, j, k, g) + wΔz(i-1, j, k, g) + wΔz(i, j-1, k, g) + wΔz(i-1, j-1, k, g); n > 0 ? s / n : Inf)
qV(i, j, k, g, f, v) = f[i, j, 1] * ℑxᶠᵃᵃ(i, j, k, g, Ayv, v) / hff(i, j, k, g)
qU(i, j, k, g, f, u) = f[i, j, 1] * ℑyᵃᶠᵃ(i, j, k, g, Axu, u) / hff(i, j, k, g)
twe_u(i, j, g, u, v) = ifelse(wetu(i, j, g), ℑyᵃᶜᵃ(i, j, 1, g, qV, fᶠᶠ, v) * Axᶠᶜᶜ(i, j, 1, g) / Vᶠᶜᶜ(i, j, 1, g), 0.0)
twe_v(i, j, g, u, v) = ifelse(wetv(i, j, g), -ℑxᶜᵃᵃ(i, j, 1, g, qU, fᶠᶠ, u) * Ayᶜᶠᶜ(i, j, 1, g) / Vᶜᶠᶜ(i, j, 1, g), 0.0)
cor = HydrostaticSphericalCoriolis(scheme = EnergyConserving())
ec_u(i, j, g, u, v) = ifelse(wetu(i, j, g), -x_f_cross_U(i, j, 1, g, cor, (u, v)), 0.0)
ec_v(i, j, g, u, v) = ifelse(wetv(i, j, g), -y_f_cross_U(i, j, 1, g, cor, (u, v)), 0.0)
u = XFaceField(grid); v = YFaceField(grid)
Random.seed!(1)
set!(u, randn(size(u)...)); set!(v, randn(size(v)...))
fill_halo_regions!((u, v))
interior(u) .*= [wetu(i, j, grid) for i in 1:size(u, 1), j in 1:size(u, 2), k in 1:1]
interior(v) .*= [wetv(i, j, grid) for i in 1:size(v, 1), j in 1:size(v, 2), k in 1:1]
fill_halo_regions!((u, v))
nu, nv = size(interior(u))[1:2], size(interior(v))[1:2]
res = Dict{Symbol, Any}()
for (name, fu, fv) in (("EnergyConserving", ec_u, ec_v), ("ThicknessWeightedEnergy", twe_u, twe_v))
    Wu = [Vᶠᶜᶜ(i, j, 1, grid) * u[i, j, 1] * fu(i, j, grid, u, v) for i in 1:nu[1], j in 1:nu[2]]
    Wv = [Vᶜᶠᶜ(i, j, 1, grid) * v[i, j, 1] * fv(i, j, grid, u, v) for i in 1:nv[1], j in 1:nv[2]]
    Wu = ifelse.(isfinite.(Wu), Wu, 0.0); Wv = ifelse.(isfinite.(Wv), Wv, 0.0)
    scale = sum(abs, Wu) + sum(abs, Wv)
    # cell map: half of each face's work to each adjacent cell
    Wc = zeros(Nx, Ny)
    for j in 1:Ny, i in 1:Nx
        Wc[i, j] = (Wu[i, j] + Wu[mod1(i + 1, Nx), j]) / 2 + (Wv[i, j] + (j < size(Wv, 2) ? Wv[i, j+1] : 0.0)) / 2
    end
    @printf("%-24s total work %+.3e, Σ|work| %.3e, ratio %.2e\n", name, sum(Wu) + sum(Wv), scale, (sum(Wu) + sum(Wv)) / scale)
    res[Symbol(name)] = Float32.(Wc)
    rowsum = vec(sum(Wu, dims = 1)) .+ [sum(Wv[:, j]) for j in 1:size(Wu, 2)]
    o = sortperm(abs.(rowsum), rev = true)[1:4]
    println("   largest row sums of work: ", join([@sprintf("j=%d: %+.2e", j, rowsum[j]) for j in o], ", "), "; rows 2:Ny-2 total ", @sprintf("%+.3e", sum(rowsum[2:end-2])))
end
jldsave(joinpath(@__DIR__, "coriolis_work.jld2"); res...)
# Corner decomposition: the work of each scheme summed by parts onto (Face, Face) corners. ThicknessWeightedEnergy:
# -f ℑx(A_y v) ℑy(A_x u)/h + f ℑy(A_x u) ℑx(A_y v)/h = 0 per corner. EnergyConserving: -f ℑx(A_y v) ℑy(Δy u) + f ℑy(A_x u) ℑx(Δx v).
Δyu(i, j, k, g, u) = ifelse(wetu(i, j, g), Δyᶠᶜᶜ(i, j, k, g) * u[i, j, k], 0.0)
Δxv(i, j, k, g, v) = ifelse(wetv(i, j, g), Δxᶜᶠᶜ(i, j, k, g) * v[i, j, k], 0.0)
corner_ec(i, j) = fᶠᶠ[i, j, 1] * (-ℑxᶠᵃᵃ(i, j, 1, grid, Ayv, v) * ℑyᵃᶠᵃ(i, j, 1, grid, Δyu, u) + ℑyᵃᶠᵃ(i, j, 1, grid, Axu, u) * ℑxᶠᵃᵃ(i, j, 1, grid, Δxv, v))
corner_tw(i, j) = (h = hff(i, j, 1, grid); isfinite(h) ? fᶠᶠ[i, j, 1] * (-ℑxᶠᵃᵃ(i, j, 1, grid, Ayv, v) * ℑyᵃᶠᵃ(i, j, 1, grid, Axu, u) + ℑyᵃᶠᵃ(i, j, 1, grid, Axu, u) * ℑxᶠᵃᵃ(i, j, 1, grid, Ayv, v)) / h : 0.0)
Cec = [corner_ec(i, j) for i in 1:Nx, j in 2:Ny]; Ctw = [corner_tw(i, j) for i in 1:Nx, j in 2:Ny]
@printf("corner maps: EnergyConserving Σ %+.3e, Σ|·| %.3e, max|·| %.3e; ThicknessWeighted Σ|·| %.3e, max|·| %.3e\n", sum(Cec), sum(abs, Cec), maximum(abs, Cec), sum(abs, Ctw), maximum(abs, Ctw))
H = [-bh[i, j, 1] for i in 1:Nx, j in 1:Ny]
jldsave(joinpath(@__DIR__, "coriolis_work_corners.jld2"); Cec = Float32.(Cec), Ctw = Float32.(Ctw), H = Float32.(H))
