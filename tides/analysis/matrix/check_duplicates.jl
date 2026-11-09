# Wet cells of the 0.36ᵒ TripolarGrid that coincide on the sphere (fold row), as the convolution operator lists them
using NumericalEarth, Oceananigans, JLD2, Printf
using Oceananigans.Grids: λnode, φnode
const NN = Base.require(Base.PkgId(Base.UUID("b8a86587-4115-5ab1-83bc-aa920d37bbce"), "NearestNeighbors"))
grid = TripolarGrid(CPU(); size = (1000, 472, 1), halo = (7, 7, 4), z = (-6000, 0))
depth = load(joinpath(@__DIR__, "..", "diag", "grid_1000.jld2"), "depth")
wet = [(i, j) for j in 1:472 for i in 1:1000 if depth[i, j] > 0]
pos(i, j) = (λ = deg2rad(λnode(i, j, 1, grid, Center(), Center(), Center())); φ = deg2rad(φnode(i, j, 1, grid, Center(), Center(), Center()));
             [cos(φ) * cos(λ), cos(φ) * sin(λ), sin(φ)])
X = reduce(hcat, [pos(i, j) for (i, j) in wet])
tree = NN.KDTree(X)
idx, dist = NN.knn(tree, X, 2, true)
d = [x[2] for x in dist]; o = [x[2] for x in idx]
for thr in (1e-12, 1e-8, 1e-6, 1e-4, 1e-3)
    @printf("pairs closer than %g rad: %d\n", thr, count(d .< thr))
end
close = findall(d .< 1e-6)
for k in close[1:min(end, 10)]
    println(wet[k], " ↔ ", wet[o[k]], " d = ", d[k], " λφ = ", round.(rad2deg.((atan(X[2, k], X[1, k]), asin(X[3, k]))), digits = 2))
end
@printf("median nearest distance %.3g rad, min %.3g; cells per j=472 row wet: %d\n", sort(d)[end ÷ 2], minimum(d), count(j == 472 for (i, j) in wet))
