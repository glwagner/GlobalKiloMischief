# Compare two dumps from dump_fields.jl: julia compare_dumps.jl a.jld2 b.jld2
using JLD2, Printf
a, b = load(ARGS[1]), load(ARGS[2])
@printf("g: %s vs %s\n", a["g"], b["g"])
for k in ("bottom", "Δzu", "Δzv", "A", "Gu0", "Gv0", "Gu4", "Gv4", "Gur", "Gvr")
    x, y = a[k], b[k]
    (isnothing(x) || isnothing(y)) && (println("$k: missing"); continue)
    d = x .- y
    i = argmax(abs.(d))
    @printf("%-6s max|a| %.4e  max|a-b| %.4e  rel %.2e  rms rel %.2e  at %s\n", k, maximum(abs, x), abs(d[i]),
            abs(d[i]) / maximum(abs, x), sqrt(sum(abs2, d) / max(sum(abs2, x), eps())), Tuple(i))
end
