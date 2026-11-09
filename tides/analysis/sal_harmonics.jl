# Harmonic analysis (M2, S2, K1, O1) of hourly η and ζ from several 1/8ᵒ runs, comparing each with the first.
# The least-squares fit is applied one snapshot at a time, so memory stays small.
#
#   julia --project=sal_env analysis/sal_harmonics.jl <first day> <last day> <reference run> <other runs...>
#
# Runs are output file names in GlobalKiloMischief/tides/output, without .jld2. Saves the complex amplitudes to
# analysis/sal_harmonics.jld2 and maps of M2 amplitude differences to analysis/sal_harmonics_M2.png.

using JLD2, LinearAlgebra, Printf, Statistics, CairoMakie

first_day, last_day = parse.(Float64, ARGS[1:2])
runs = ARGS[3:end]
output = joinpath(@__DIR__, "..", "GlobalKiloMischief", "tides", "output")

periods = (M2 = 12.4206012, S2 = 12.0, K1 = 23.9344696, O1 = 25.8193417) # hours
ω = [2π / (3600T) for T in values(periods)]
H = 7 # halo

"""Complex amplitudes Z (one array per constituent) with field ≈ c + Σ Re(Zₖ exp(iωₖt))."""
function harmonic_fit(run, name)
    jldopen(joinpath(output, run * ".jld2")) do file
        iterations = keys(file["timeseries/t"])
        times = [file["timeseries/t/$i"] for i in iterations]
        window = findall(t -> first_day * 86400 ≤ t ≤ last_day * 86400, times)
        τ = times[window]
        A = hcat(ones(length(τ)), [cos(ωₖ * t) for t in τ, ωₖ in ω], [sin(ωₖ * t) for t in τ, ωₖ in ω])
        P = pinv(A) # (1 + 2K) × length(τ)
        K = length(ω)
        coefficients = nothing
        for (n, i) in enumerate(iterations[window])
            snapshot = Float64.(file["timeseries/$name/$i"][H+1:end-H, H+1:end-H, 1])
            isnothing(coefficients) && (coefficients = [zeros(size(snapshot)) for _ in 1:2K+1])
            for q in 1:2K+1
                coefficients[q] .+= P[q, n] .* snapshot
            end
        end
        Z = [coefficients[1 + k] .- im .* coefficients[1 + K + k] for k in 1:K]
        return NamedTuple{keys(periods)}(Z), length(window)
    end
end

amplitudes = Dict{String, Any}()
for run in runs, name in ("η", "ζ")
    t = @elapsed Z, nt = harmonic_fit(run, name)
    amplitudes["$run/$name"] = Z
    @info @sprintf("%s %s: fitted %d snapshots in %.0f s", run, name, nt, t)
end

reference = first(runs)
# ζ is constant (-650 m) where the ocean is shallower than its level, so each field gets its own mask
wets = Dict(name => abs.(amplitudes["$reference/$name"].M2) .> 0 for name in ("η", "ζ"))
wet = wets["η"]
Ny = size(wet, 2)
φ = [-80 + 170 * (j - 0.5) / Ny for i in 1:size(wet, 1), j in 1:Ny] # exact south of the fold
area = cosd.(clamp.(φ, -90, 90)) # approximate area weights; the fold region is weighted as if regular
area_rms(x, wet) = sqrt(sum(area[wet] .* abs2.(x[wet])) / sum(area[wet]))

for name in ("η", "ζ"), c in keys(periods)
    local wet = wets[name]
    rms(x) = area_rms(x, wet)
    Zʳ = amplitudes["$reference/$name"][c]
    @printf("%-58s %s %s: rms amplitude %.4f m\n", reference, name, c, rms(Zʳ))
    for run in runs[2:end]
        Z = amplitudes["$run/$name"][c]
        @printf("%-58s %s %s: rms amplitude %.4f m, rms |ΔZ| %.4f m (%.2f of reference), max |ΔZ| %.2f m\n",
                run, name, c, rms(Z), rms(Z .- Zʳ), rms(Z .- Zʳ) / rms(Zʳ), maximum(abs, (Z .- Zʳ)[wet]))
    end
end

jldsave(joinpath(@__DIR__, "sal_harmonics.jld2"); runs, first_day, last_day,
        amplitudes = Dict(k => Dict(String(c) => v[c] for c in keys(periods)) for (k, v) in amplitudes))

fig = Figure(size = (1700, 400 * length(runs)))
for (r, run) in enumerate(runs), (col, name) in enumerate(("η", "ζ"))
    local wet = wets[name]
    masked(x) = ifelse.(wet, x, NaN)
    Z = amplitudes["$run/$name"].M2
    if r == 1
        field, title, colormap = abs.(Z), "M2 |$name| (m), $run", :viridis
        limit = quantile(field[wet], 0.99)
        range = (0, limit)
    else
        field, title, colormap = abs.(Z) .- abs.(amplitudes["$reference/$name"].M2), "M2 |$name| − reference (m), $run", :balance
        limit = quantile(abs.(field[wet]), 0.99)
        range = (-limit, limit)
    end
    ax = Axis(fig[r, 2col - 1]; title, titlesize = 11, aspect = DataAspect())
    hm = heatmap!(ax, masked(field); colormap, colorrange = range)
    Colorbar(fig[r, 2col], hm)
    hidedecorations!(ax)
end
save(joinpath(@__DIR__, "sal_harmonics_M2.png"), fig)
@info "Saved analysis/sal_harmonics_M2.png and analysis/sal_harmonics.jld2"
