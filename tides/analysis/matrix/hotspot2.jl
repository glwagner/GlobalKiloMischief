using JLD2, Printf
depth, λ, φ = load(joinpath(@__DIR__, "..", "diag", "grid_1000.jld2"), "depth", "λ", "φ")
out = joinpath(@__DIR__, "..", "..", "GlobalKiloMischief", "tides", "output")
h = load(joinpath(out, ARGS[1] * "_harmonics.jld2"))
m = h["mean"]
@printf("(640,50): mean η %.2f m; constituents: %s\n", m[640, 50], join([@sprintf("%s %.2f", k, abs(v[640, 50])) for (k, v) in h["Z"]], ", "))
for k in sort(collect(keys(h["η_snapshots"]))); s = h["η_snapshots"][k]; @printf("  %s η(640,50) = %.2f, max|η - mean| over cells %.2f\n", k, s[640, 50], maximum(abs, filter(isfinite, s .- m))); end
I, J = 712:728, 446:462
for (lab, F) in (("depth", depth), ("|M2|", abs.(h["Z"]["M2"])), ("λ", λ), ("φ", φ))
    println(lab); for j in reverse(J); println(j, ": ", join([@sprintf("%5.0f", lab == "|M2|" ? 100F[i, j] : F[i, j]) for i in I], " ")); end
end
println("mean η around (640,50)")
for j in 54:-1:46; println(j, ": ", join([@sprintf("%6.2f", m[i, j]) for i in 636:644], " ")); end
o = sortperm(vec(abs.(m)), rev = true)[1:15]
for q in o; i, j = Tuple(CartesianIndices(m)[q]); @printf("  mean %.2f m at (%d,%d) %.1fE %.1fN H %.0f, |M2| %.2f\n", m[q], i, j, λ[q], φ[q], depth[q], abs(h["Z"]["M2"][q])); end
