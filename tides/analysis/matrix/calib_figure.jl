# (τ, χ) error surface of the calibration sweep (EOT20 deep, plus held-out TPXO9 deep and gauges > 200 m)
using CairoMakie, DelimitedFiles, Printf
d = readdlm(joinpath(@__DIR__, "calib_scores.tsv"), '\t'; header = true)[1]   # χ, τ_days, EOT_deep, TPXO_deep, gauge_200
χ, τ, e, t, g = Float64.(d[:, 1]), Float64.(d[:, 2]), Float64.(d[:, 3]), Float64.(d[:, 4]), Float64.(d[:, 5])
fig = Figure(size = (1500, 480), fontsize = 14)
for (k, (v, lab)) in enumerate(((e, "EOT20 deep (tuning)"), (t, "TPXO9 deep (held out)"), (g, "TICON-3 gauges > 200 m (held out)")))
    ax = Axis(fig[1, 2k-1]; xscale = log10, xlabel = "biharmonic timescale τ (days)", ylabel = "wave-drag scale χ (× JSL_RH)", title = lab,
              xticks = [0.5, 1, 3, 10, 30, 100])
    sc = scatter!(ax, τ, χ; color = v, colormap = :viridis, markersize = 34, marker = :rect)
    for i in eachindex(v); text!(ax, τ[i], χ[i]; text = @sprintf("%.2f", v[i]), align = (:center, :center), fontsize = 10, color = :white); end
    i = argmin(v); scatter!(ax, [τ[i]], [χ[i]]; color = :transparent, strokecolor = :gray50, strokewidth = 2, markersize = 40, marker = :rect)
    scatter!(ax, [1.0], [1.0]; color = :transparent, strokecolor = :red, strokewidth = 3, markersize = 44, marker = :rect)
    Colorbar(fig[1, 2k], sc; label = "rms |ΔZ|/√2 (cm)")
end
Label(fig[0, :], "Calibration sweep: MOM6 aw depth (min 5 m), JSL_RH × χ, flat biharmonic A = Δ⁴/τ, thickness-weighted Coriolis, Δt 5 min, SH n=40 PREM, M2 only\nred: chosen χ = 1.0, τ = 1 d (χ tuned on EOT20 deep only; τ chosen for stability margin (4.5× the bound at Δt 5 min) and conservatism, not on the held-out gauges); gray: minimum of each panel")
save(joinpath(@__DIR__, "calibration_surface.png"), fig)
