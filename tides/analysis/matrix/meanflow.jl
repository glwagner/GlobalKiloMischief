# Mean (fitted) η and velocities: where are they large?
using JLD2, Printf, Statistics
depth, λ, φ = load(joinpath(@__DIR__, "..", "diag", "grid_1000.jld2"), "depth", "λ", "φ")
out = joinpath(@__DIR__, "..", "..", "GlobalKiloMischief", "tides", "output")
for name in ARGS
    h = load(joinpath(out, name * "_harmonics.jld2")); m = h["mean"]
    @printf("== %s: max|mean η| %.2f m, area-unweighted rms %.3f m, cells |mean η|>1 m: %d\n", name, maximum(abs, m), sqrt(mean(abs2, m[depth .> 0])), count(abs.(m) .> 1))
    if haskey(h, "mean_u")
        u, v = h["mean_u"], h["mean_v"]; s = sqrt.(u[1:size(v,1), 1:size(u,2)] .^ 2 .+ v[1:size(v,1), 1:size(u,2)] .^ 2)
        o = sortperm(vec(s), rev = true)[1:8]
        for q in o; i, j = Tuple(CartesianIndices(s)[q]); @printf("  mean speed %.2f m/s (u %.2f v %.2f) at (%d,%d) %.1fE %.1fN H %.0f, mean η %.2f, |M2 u| %.2f\n", s[q], u[i,j], v[i,j], i, j, λ[i,j], φ[i,j], depth[i,j], m[i,j], abs(h["Zu"]["M2"][i,j])); end
        @printf("  cells with mean speed > 0.5 m/s: %d, > 0.1: %d\n", count(s .> 0.5), count(s .> 0.1))
    end
end
