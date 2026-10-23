# Where does η grow before the blow-up, and how deep is the ocean there?
using Oceananigans
using Oceananigans.Grids: λnodes, φnodes
η = FieldTimeSeries(joinpath(@__DIR__, "output", ARGS[1] * ".jld2"), "η"; backend = InMemory())
grid = η.grid
Nx, Ny = size(grid, 1), size(grid, 2)
zb = Array(parent(grid.immersed_boundary.bottom_height))[1 .+ grid.Hx:Nx + grid.Hx, 1 .+ grid.Hy:Ny + grid.Hy, 1]
λ = λnodes(grid, Center(), Center(), Center())[1:Nx, 1:Ny]
φ = φnodes(grid, Center(), Center(), Center())[1:Nx, 1:Ny]
for n in eachindex(η.times)
    ηⁿ = Array(interior(η[n], :, :, 1))
    any(isnan, ηⁿ) && (println("NaN at hour ", η.times[n] / 3600); break)
    n % 6 == 0 || continue
    value, index = findmax(abs, ηⁿ)
    i, j = Tuple(index)
    println("hour ", round(Int, η.times[n] / 3600), ": max|η| = ", round(value, digits=2), " m at (i, j) = ", (i, j),
            ", (λ, φ) = ", round.((λ[i, j], φ[i, j]), digits=1), ", depth ", round(-zb[i, j]), " m; min(depth + η) = ",
            round(minimum((ηⁿ .- zb)[zb .< 0]), digits=2), " m")
end
