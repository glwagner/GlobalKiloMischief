# 2×4 matrix figure (EOT20 M2 amplitude alone on top; panels |Z − Z_EOT20| with scores in titles) and bar charts of
# errors and cost, from matrix/runs_latest.tsv, scores_latest.tsv, costs_latest.tsv. Writes matrix_amplitude.png,
# matrix_bars.png in matrix/.
using JLD2, CairoMakie, Printf, Dates, DelimitedFiles
const NCDatasets = Base.require(Base.PkgId(Base.UUID("85f8d34a-cbdd-5861-8df4-14fed0d494ab"), "NCDatasets"))
here = @__DIR__
output = joinpath(here, "..", "..", "GlobalKiloMischief", "tides", "output")
depth, λ, φ = load(joinpath(here, "..", "diag", "grid_1000.jld2"), "depth", "λ", "φ")
mom6depth = load(joinpath(here, "..", "..", "GlobalKiloMischief", "tides", "mom6_p04arc9_on_1000.jld2"), "depth")
runs = readdlm(joinpath(here, "runs_latest.tsv"), '\t', String; comments = true)        # row, column, run name, note
sc = readdlm(joinpath(here, "scores_latest.tsv"), '\t'; header = true)
num(x) = x isa Number ? Float64(x) : NaN
scores = Dict(String(r[1]) => [num(x) for x in r[2:end]] for r in eachrow(sc[1]))
col(name) = findfirst(==(name), vec(sc[2])) - 1
co = readdlm(joinpath(here, "costs_latest.tsv"), '\t'; header = true)
costs = Dict(String(r[1]) => [num(x) for x in r[2:end]] for r in eachrow(co[1]))
rows = unique(runs[:, 1]); cols = unique(runs[:, 2])
data = joinpath(here, "..", "..", "egbert_ray", "data")
eot = NCDatasets.NCDataset(joinpath(data, "EOT20/ocean_tides/M2_ocean_eot20.nc"))
lo, la = Float64.(eot["lon"][:]), Float64.(eot["lat"][:])
EZ = map((a, g) -> ismissing(a) || ismissing(g) ? complex(NaN) : a / 100 * cis(-deg2rad(g)), eot["amplitude"][:, :], eot["phase"][:, :]); close(eot)
function bilinear(x0, y0)
    x = mod(x0 - lo[1], 360) / (lo[2] - lo[1]) + 1; y = (y0 - la[1]) / (la[2] - la[1]) + 1
    i = clamp(floor(Int, x), 1, length(lo) - 1); j = clamp(floor(Int, y), 1, length(la) - 1); a, b = x - i, y - j
    (1 - a) * (1 - b) * EZ[i, j] + a * (1 - b) * EZ[i+1, j] + (1 - a) * b * EZ[i, j+1] + a * b * EZ[i+1, j+1]
end
EOT = bilinear.(λ, φ)
jc(d) = (Dates.datetime2julian(d) - 2451545.0) / 36525
T₀ = jc(DateTime(2015, 1, 1)); V₀ = mod(2(280.4664567 + 36000.76983T₀) - 2(218.3164477 + 481267.88123421T₀), 360)
function greenwich(Z, mid)
    N = 125.04452 - 1934.136261jc(DateTime(2015, 1, 1) + Millisecond(round(Int, mid * 86400e3)))
    u = -2.14sind(N); f = 1.0004 - 0.0373cosd(N) + 0.0002cosd(2N)
    @. abs(Z) / f * cis(angle(Z) - deg2rad(V₀ + u))
end
model(name) = (h = load(joinpath(output, name * "_harmonics.jld2")); greenwich(h["Z"]["M2"], (h["first_day"] + h["last_day"]) / 2))
wet = (depth .> 0) .& (mom6depth .> 0) .& isfinite.(EOT)
fig = Figure(size = (2400, 1500), fontsize = 14)
nc = length(cols)
ax = Axis(fig[1, 2:3]; title = "EOT20 M2 amplitude (truth)", aspect = DataAspect()); hidedecorations!(ax)
hm = heatmap!(ax, ifelse.(wet, abs.(EOT), NaN); colormap = :viridis, colorrange = (0, 1.5), nan_color = :gray85)
Colorbar(fig[1, nc + 1], hm; label = "M2 amplitude (m)", height = Relative(0.8))
for (r, rowname) in enumerate(rows), (c, colname) in enumerate(cols)
    k = findfirst(i -> runs[i, 1] == rowname && runs[i, 2] == colname, 1:size(runs, 1)); isnothing(k) && continue
    name, note = runs[k, 3], runs[k, 4]
    s = get(scores, name, fill(NaN, 20))
    ax = Axis(fig[r + 1, c]; aspect = DataAspect(), titlesize = 13,
              title = @sprintf("%s, %s%s\nEOT20 deep %.2f (amp %.2f) · TPXO9 deep %.2f · gauges %.1f / %.1f cm", rowname, colname,
                               isempty(note) ? "" : " ($note)", s[col("EOT_deep")], s[col("EOT_deep_amp")], s[col("TPXO_deep")], s[col("gauge_all")], s[col("gauge_200")]))
    hidedecorations!(ax)
    if isfile(joinpath(output, name * "_harmonics.jld2"))
        E = abs.(model(name) .- EOT)
        global eh = heatmap!(ax, ifelse.(wet, E, NaN); colormap = :inferno, colorrange = (0, 0.3), highclip = :white, nan_color = :gray85)
    else
        heatmap!(ax, ifelse.(wet, 0.0, NaN); colormap = [:gray97, :gray97], nan_color = :gray85)
        text!(ax, size(depth, 1) / 2, size(depth, 2) / 2; text = "pending", align = (:center, :center), fontsize = 30, color = :gray40)
        ax.title = "$rowname, $colname ($note)"
    end
end
Colorbar(fig[2:3, nc + 1], eh; label = "|Z − Z_EOT20| (m)")
Label(fig[0, :], "M2 complex error vs EOT20, single-layer 0.36ᵒ model, \"calibrated\" numerics (χ 1.0 JSL_RH, biharmonic τ 1 d, thickness-weighted Coriolis, min depth 5 m, Δt 5 min, M2 only, fit days 30–60); errors rms|ΔZ|/√2 (cm): deep = H > 1000 m, |φ| < 66ᵒ; gauges all / > 200 m",
      fontsize = 16)
save(joinpath(here, "matrix_amplitude.png"), fig; px_per_unit = 1)

fig = Figure(size = (1500, 900), fontsize = 13)
labels = [@sprintf("%s\n%s", runs[i, 1], runs[i, 2]) for i in 1:size(runs, 1)]
names = runs[:, 3]
metrics = [("EOT_deep", "EOT20 deep"), ("TPXO_deep", "TPXO9 deep"), ("EOT_global", "EOT20 global"), ("gauge_200", "gauges > 200 m")]
ax = Axis(fig[1, 1]; limits = (nothing, (0, 1.25maximum(filter(isfinite, [get(scores, n, fill(NaN, 20))[col("gauge_200")] for n in names])))),
          ylabel = "rms |ΔZ|/√2 (cm)", title = "M2 errors", xticks = (1:length(names), labels), xticklabelsize = 10)
colors = Makie.wong_colors()
for (q, (m, l)) in enumerate(metrics)
    barplot!(ax, (1:length(names)) .+ (q - 2.5) * 0.2, [get(scores, n, fill(NaN, 20))[col(m)] for n in names]; width = 0.19, color = colors[q], label = l)
end
axislegend(ax; position = :rt, nbanks = 4)
ax2 = Axis(fig[2, 1]; ylabel = "wall s per model day (T4)", yscale = log10, title = "Cost on one T4: whole run and SAL overhead vs the no-SAL run of the same row (≤ 0.1 shown at 0.1)",
           xticks = (1:length(names), labels), xticklabelsize = 10)
tot = [get(costs, n, fill(NaN, 3))[1] for n in names]
noneof(i) = names[findfirst(k -> runs[k, 1] == runs[i, 1] && runs[k, 2] == "none", 1:length(names))]
sal = [tot[i] - get(costs, noneof(i), fill(NaN, 3))[1] for i in eachindex(names)]
barplot!(ax2, (1:length(names)) .- 0.2, tot; fillto = 0.1, width = 0.38, color = colors[1], label = "whole run")
barplot!(ax2, (1:length(names)) .+ 0.2, max.(ifelse.(isfinite.(sal), sal, 0.1), 0.1); fillto = 0.1, width = 0.38, color = colors[2], label = "SAL overhead")
axislegend(ax2; position = :rt)
save(joinpath(here, "matrix_bars.png"), fig; px_per_unit = 1.5)
