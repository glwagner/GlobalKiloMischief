# Frequency content of hourly η from the single-layer runs: area-weighted power in tidal and non-tidal bands over a
# window, and the global mean η² every day (to detect spurious growing modes).
#   julia --project=sal_env analysis/barotropic_spectra.jl <run suffix>... (e.g. none scalar spherical_harmonic)
using JLD2, Printf, Statistics
const FFTW = Base.require(Base.PkgId(Base.UUID("7a1cc6ca-52ef-59f5-83cd-3a7055c09341"), "FFTW")); rfft(x, d) = FFTW.rfft(x, d)
using Oceananigans, NumericalEarth
using Oceananigans.Operators: Azᶜᶜᶜ
output = joinpath(@__DIR__, "..", "GlobalKiloMischief", "tides", "output")
grid = TripolarGrid(CPU(); size = (1000, 472, 1), halo = (7, 7, 4), z = (-6000, 0))
area = [Azᶜᶜᶜ(i, j, 1, grid) for i in 1:1000, j in 1:472]
for sal in ARGS
    jldopen(joinpath(output, "barotropic_tides_0.36deg_dt10min_$sal.jld2")) do file
        its = keys(file["timeseries/t"]); t = [file["timeseries/t/$i"] for i in its] ./ 86400
        η = [Float64.(file["timeseries/η/$i"][8:1007, 8:479, 1]) for i in its]
        wet = η[end] .!= 0
        w = area[wet] ./ sum(area[wet])
        daily = [@sprintf("%.3f", sqrt(sum(w .* η[n][wet] .^ 2))) for n in eachindex(t) if isinteger(t[n]) && t[n] > 0]
        window = findall(x -> 5 ≤ x < 20, t)
        X = reduce(hcat, [η[n][wet] for n in window])          # points × hours
        X .-= mean(X, dims = 2)
        P = abs2.(rfft(X, 2)) ./ size(X, 2)^2 .* 2
        f = (0:size(P, 2)-1) ./ (size(X, 2) / 24)                # cycles per day
        band(lo, hi) = sqrt(sum(w .* sum(P[:, lo .≤ f .< hi], dims = 2)))
        @info @sprintf("%-20s rms η by band, days 5–20: diurnal %.3f m, semidiurnal %.3f m, 2.5–3.5 cpd %.3f m, 3.5–4.5 cpd %.3f m, 4.5–12 cpd %.3f m, <0.7 cpd %.3f m",
                       sal, band(0.7, 1.3), band(1.7, 2.3), band(2.5, 3.5), band(3.5, 4.5), band(4.5, 12.1), band(0.01, 0.7))
        @info "$sal: global rms η at the end of each day (m): " * join(daily, " ")
    end
end
