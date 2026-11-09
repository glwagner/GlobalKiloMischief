using JLD2, Printf, Statistics
H = 7
jldopen(ARGS[1], "r") do f
    its = sort(parse.(Int, keys(f["timeseries/t"])))
    λ = [(i - 0.5) / 8 for i in 1:2880, j in 1:1360]; φ = [-80 + (j - 0.5) / 8 for i in 1:2880, j in 1:1360]  # approximate (offset in λ unknown; exact south of the fold)
    for i in its[end-8:4:end]
        t = f["timeseries/t/$i"] / 86400
        for v in ("ζ", "w", "η")
            x = f["timeseries/$v/$i"][H+1:end-H, H+1:end-H, 1]
            wet = v == "ζ" ? abs.(x .+ 650) .> 1 : trues(size(x))
            a = abs.(x) .* wet
            m, I = findmax(a)
            q = quantile(a[wet], [0.5, 0.99, 0.999])
            @printf("day %.2f %s: max %.3g at (i,j)=(%d,%d) λ=%.1f φ=%.1f; median %.2g p99 %.2g p99.9 %.2g\n", t, v, m, I[1], I[2], λ[I], φ[I], q...)
        end
    end
end
