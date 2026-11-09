# Time ConvolutionLoadingOperator(accurate) on the immersed grids of the two configurations
using NumericalEarth, Oceananigans, CUDA, JLD2, Printf
using Oceananigans.BoundaryConditions: fill_halo_regions!
const O = NumericalEarth.Oceans
depth = load(joinpath(@__DIR__, "..", "diag", "grid_1000.jld2"), "depth")
m6 = load(joinpath(@__DIR__, "..", "..", "GlobalKiloMischief", "tides", "mom6_p04arc9_on_1000.jld2"), "depth")
for (label, D, Hmax, minimum_depth) in (("ours", depth, 6000.0, 20.0), ("mom6", m6, 9314.0, 5.0))
    ug = TripolarGrid(GPU(); size = (1000, 472, 1), halo = (7, 7, 4), z = (-Hmax, 0))
    zb = Field{Center, Center, Nothing}(ug); set!(zb, reshape(map(x -> x > 0 ? -max(x, minimum_depth) : 0.0, D), 1000, 472, 1)); fill_halo_regions!(zb)
    grid = ImmersedBoundaryGrid(ug, PartialCellBottom(zb; minimum_fractional_cell_height = min(1e-3, minimum_depth / 2Hmax)); active_cells_map = true)
    for kw in ((; FT = Float32), (; FT = Float32, remove_duplicates = true, subtract_uniform = true))
        op = O.ConvolutionLoadingOperator(grid; components = :grid, kw...)
        p = Field{Center, Center, Nothing}(grid); set!(p, (λ, φ) -> cosd(φ)^2 * cosd(2λ))
        a, b = Field{Center, Center, Nothing}(grid), Field{Center, Center, Nothing}(grid)
        O.loading_gradient!(a, b, op, p); CUDA.synchronize()
        t = @elapsed (O.loading_gradient!(a, b, op, p); CUDA.synchronize())
        @printf("%s %s: N = %d, %.2f s per application, any NaN %s\n", label, kw, length(op.x), t, any(isnan, Array(interior(a))))
        flush(stdout)
    end
end
