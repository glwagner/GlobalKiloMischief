# Extra figures: Coriolis energy test, Δt convergence vs the explicit adjoint-session model, mean η near Antarctica
# before/after the thickness-weighted Coriolis.
using CairoMakie, JLD2, Printf
here = @__DIR__; logs = joinpath(here, "..", "..", "logs")
out = joinpath(here, "..", "..", "GlobalKiloMischief", "tides", "output")
function energy(file)
    d = Float64[]; E = Float64[]
    for l in eachline(joinpath(logs, file))
        m = match(r"ENERGY day ([0-9.]+): K \S+ P \S+ E ([0-9.e+]+)", l); m === nothing && continue
        push!(d, parse(Float64, m[1])); push!(E, parse(Float64, m[2]))
    end
    d, E
end
fig = Figure(size = (900, 520), fontsize = 15)
ax = Axis(fig[1, 1]; yscale = log10, xlabel = "day", ylabel = "total energy E / E(day 0.25)",
          title = "Unforced, inviscid, no drag; native depth, 5 m minimum; seeded with 0.05 m/s noise; Δt 10 min")
for (f, lab, c) in (("etEnergyConserving-3988.log", "Oceananigans EnergyConserving (blows up at day ≈ 4.5)", :orangered3),
                    ("etThicknessWeightedEnergy-3989.log", "ThicknessWeightedEnergy (corner PV f/h, face h)", :dodgerblue4))
    d, E = energy(f); k = findfirst(>(0.2), d)
    scatterlines!(ax, d[k:end], E[k:end] ./ E[k]; color = c, label = lab)
end
axislegend(ax; position = :lt, framevisible = false)
save(joinpath(here, "coriolis_energy_test.png"), fig)

# Δt convergence (old Coriolis on both sides; no-SAL, τ 30 d, χ 0.8; deep H > 1000 m, |φ| < 66ᵒ; from dtcmp-3991.log)
Δt = [10, 5, 2.5]; ratio = 1 ./ [1.0279, 1.0140, 1.0072]; rmsd = [1.17, 0.65, 0.43]; dph = [1.22, 0.77, 0.62]
fig = Figure(size = (1100, 450), fontsize = 15)
a = Axis(fig[1, 1]; xlabel = "Δt (min)", ylabel = "deep mean amplitude, split-explicit / explicit", title = "M2 amplitude vs time step",
         xticks = [0, 2.5, 5, 10])
scatterlines!(a, Δt, ratio; color = :dodgerblue4, markersize = 12); scatter!(a, [0], [1]; color = :black, marker = :star5, markersize = 18)
lines!(a, [0, 10], [1, ratio[1]]; color = :gray, linestyle = :dash); xlims!(a, -0.5, 10.5)
text!(a, 0.3, 1.0; text = "explicit, Δt 22.9 s", align = (:left, :top))
b = Axis(fig[1, 2]; xlabel = "Δt (min)", ylabel = "deep rms|ΔZ|/√2 (cm) / mean |Δphase| (ᵒ)", title = "Difference from the explicit model", xticks = [2.5, 5, 10])
scatterlines!(b, Δt, rmsd; label = "rms|ΔZ|/√2 (cm)", color = :dodgerblue4); scatterlines!(b, Δt, dph; label = "|Δphase| (ᵒ)", color = :orangered3, marker = :rect)
axislegend(b; position = :lt, framevisible = false)
Label(fig[0, :], "Split-explicit free-surface damping is linear in Δt (averaging kernel; substep count has no effect). Calibration uses Δt 5 min.")
save(joinpath(here, "dt_convergence.png"), fig)

# mean η near Antarctica: before (default EnstrophyConserving Coriolis) and after (ThicknessWeightedEnergy)
λ, φ, depth = load(joinpath(here, "..", "diag", "grid_1000.jld2"), "λ", "φ", "depth")
J = findall(j -> φ[1, j] < -50, 1:size(φ, 2))
I = 330:700
fig = Figure(size = (1700, 560), fontsize = 15)
for (k, (n, lab)) in enumerate((("barotropic_tides_0.36deg_dt10min_none_wavedrag1.5_s0_ref", "before: Oceananigans Coriolis (EnstrophyConserving)"),
                                ("barotropic_tides_0.36deg_dt10min_none_wavedrag1.5_ThicknessWeightedEnergy_s0_twe", "after: ThicknessWeightedEnergy")))
    m = load(joinpath(out, n * "_harmonics.jld2"), "mean")
    ax = Axis(fig[1, k]; title = @sprintf("%s\nmax |mean η| %.2f m", lab, maximum(abs, m[:, J])), aspect = AxisAspect(1.6))
    hidedecorations!(ax)
    global hm = heatmap!(ax, I, J, ifelse.(depth[I, J] .> 0, m[I, J], NaN); colormap = :balance, colorrange = (-3, 3), nan_color = :gray85)
end
Colorbar(fig[1, 3], hm; label = "fitted mean η, days 20–30 (m)")
Label(fig[0, :], "Mean sea level south of 50ᵒS, grid columns i = 330–700 (≈ 190ᵒE–320ᵒE: Ross Sea to the Antarctic Peninsula; native depth, 5 m minimum, no SAL, M2 only, 30 days)")
save(joinpath(here, "antarctic_mean_eta.png"), fig)
