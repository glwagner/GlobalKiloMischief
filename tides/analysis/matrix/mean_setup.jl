# Is the steady η (fitted mean) the set-up η' = -(b/g) H/2 + c expected from a uniform buoyancy b = g(αT - βS) in
# Oceananigans' hydrostatic pressure at partial-cell centres (single layer, cell centre at -H/2)?
using JLD2, Printf, Statistics
depth = load(joinpath(@__DIR__, "..", "diag", "grid_1000.jld2"), "depth")
out = joinpath(@__DIR__, "..", "..", "GlobalKiloMischief", "tides", "output")
for name in ARGS
    h = load(joinpath(out, name * "_harmonics.jld2"))
    m = h["mean"]; H = haskey(ENV, "MOM6") ? load(joinpath(out, "..", "mom6_p04arc9_on_1000.jld2"), "depth") : depth
    w = (H .> 0) .& isfinite.(m) .& (m .!= 0)
    x, y = H[w], m[w]
    slope = cov(x, y) / var(x); r = cor(x, y)
    @printf("%s: mean η range %.2f to %.2f m, std %.2f; regression on H: slope %.5f (r = %.3f)\n", name, extrema(y)..., std(y), slope, r)
    @printf("   predicted slope -(b/g)/2 with b = g(1.67e-4·10 - 7.8e-4·35) = %.5f\n", -(1.67e-4 * 10 - 7.8e-4 * 35) / 2)
end
