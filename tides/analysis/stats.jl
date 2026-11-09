# Daily max|η|, max|ζ|, max|w| and NaN count from a JLD2Writer output file (no Oceananigans needed)
using JLD2, Printf
jldopen(ARGS[1], "r") do f
    its = sort(parse.(Int, keys(f["timeseries/t"])))
    println(length(its), " snapshots, last t = ", f["timeseries/t/$(last(its))"] / 86400, " days")
    println("size ζ: ", size(f["timeseries/ζ/$(last(its))"]), " eltype ", eltype(f["timeseries/ζ/$(last(its))"]))
    for i in its
        t = f["timeseries/t/$i"] / 86400
        (isinteger(round(t * 24)) && round(Int, t * 24) % parse(Int, get(ENV, "EVERY", "24")) == 0) || i == last(its) || continue
        s = Any[f["timeseries/$v/$i"] for v in ("η", "ζ", "w")]
        s[2] = s[2][abs.(s[2] .+ 650) .> 1]  # immersed cells have T = 0, so ζ = -d/2 = -650 m there
        nans = sum(count(isnan, x) for x in s)
        @printf("day %5.2f  max|η| = %.2f m  max|ζ| = %.1f m  max|w| = %.1e m/s  NaN = %d\n", t, (maximum(abs, filter(!isnan, x)) for x in s)..., nans)
    end
end
