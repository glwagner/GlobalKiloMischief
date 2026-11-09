using JLD2, Printf, Statistics
d = load(joinpath(@__DIR__, "..", "m2_gauges.jld2"))
println(keys(d["at_gauges"])); println(keys(d["gauges"]))
for (gk, g) in d["gauges"], (rk, r) in d["at_gauges"]
    haskey(r, gk) || continue
    m = abs.(r[gk]); o = abs.(g.Z)
    rk == "0.36° 1L SH n=40 χ=3" || continue
    sel = findall(m .> 1.8)
    println("== $gk / $rk: $(length(sel)) of $(length(m)) gauges with model > 1.8 m")
    for i in sel[sortperm(m[sel], rev=true)][1:min(end, 40)]
        @printf("  φ=%6.2f λ=%7.2f model %.2f m obs %.2f m depth %.0f\n", g.φ[i], g.λ[i], m[i], o[i], get(d["depth_at"], gk, fill(NaN, length(m)))[i])
    end
end
