# Figures for matrix/RESULTS.md: EOT20 M2 amplitude beside 2×4 maps of the complex M2 error |Z_model − Z_EOT20| (rows:
# configurations, columns: SAL methods) with errors in the titles, and bar charts of the errors and of the cost.
#   julia --project=sal_env analysis/matrix/figures.jl   (reads matrix/runs.tsv, matrix/scores.tsv, matrix/costs.tsv)
using JLD2, CairoMakie, Printf, Dates, DelimitedFiles
const NCDatasets = Base.require(Base.PkgId(Base.UUID("85f8d34a-cbdd-5861-8df4-14fed0d494ab"), "NCDatasets"))
here = @__DIR__
output = joinpath(here, "..", "..", "GlobalKiloMischief", "tides", "output")
depth, λ, φ = load(joinpath(here, "..", "diag", "grid_1000.jld2"), "depth", "λ", "φ")
mom6depth = load(joinpath(here, "..", "..", "GlobalKiloMischief", "tides", "mom6_p04arc9_on_1000.jld2"), "depth")
runs = readdlm(joinpath(here, "runs.tsv"), '\t', String; comments = true)       # row, column, run name
sc = readdlm(joinpath(here, "scores.tsv"), '\t'; header = true)
scores = Dict(String(r[1]) => [x isa Number ? Float64(x) : NaN for x in r[2:end]] for r in eachrow(sc[1]))
col(name) = findfirst(==(name), vec(sc[2])) - 1
co = readdlm(joinpath(here, "costs.tsv"), '\t'; header = true)
costs = Dict(String(r[1]) => [x isa Number ? Float64(x) : NaN for x in r[2:end]] for r in eachrow(co[1]))   # s per day, ms per step, SAL ms per update, steps per day
rows = unique(runs[:, 1]); cols = unique(runs[:, 2])
jc(d) = (Dates.datetime2julian(d) - 2451545.0) / 36525
# EOT20 at cell centres and model M2 in Greenwich phase, as in score_matrix.jl
eot = NCDatasets.NCDataset(joinpath(here, "..", "..", "egbert_ray", "data", "EOT20", "ocean_tides", "M2_ocean_eot20.nc"))
lo, la = Float64.(eot["lon"][:]), Float64.(eot["lat"][:])
EZ = map((a, g) -> ismissing(a) || ismissing(g) ? complex(NaN) : a / 100 * cis(-deg2rad(g)), eot["amplitude"][:, :], eot["phase"][:, :]); close(eot)
function bilinear(x0, y0)
    x = mod(x0 - lo[1], 360) / (lo[2] - lo[1]) + 1; y = (y0 - la[1]) / (la[2] - la[1]) + 1
    i = clamp(floor(Int, x), 1, length(lo) - 1); j = clamp(floor(Int, y), 1, length(la) - 1); a, b = x - i, y - j
    (1 - a) * (1 - b) * EZ[i, j] + a * (1 - b) * EZ[i+1, j] + (1 - a) * b * EZ[i, j+1] + a * b * EZ[i+1, j+1]
end
EOT = bilinear.(λ, φ)
T₀ = jc(DateTime(2015, 1, 1)); V₀ = mod(2(280.4664567 + 36000.76983T₀) - 2(218.3164477 + 481267.88123421T₀), 360)
function greenwich(Z, mid)
    N = 125.04452 - 1934.136261jc(DateTime(2015, 1, 1) + Millisecond(round(Int, mid * 86400e3)))
    u = -2.14sind(N); f = 1.0004 - 0.0373cosd(N) + 0.0002cosd(2N)
    @. abs(Z) / f * cis(angle(Z) - deg2rad(V₀ + u))
end
error_magnitude(name) = (h = load(joinpath(output, name * "_harmonics.jld2")); abs.(greenwich(h["Z"]["M2"], (h["first_day"] + h["last_day"]) / 2) .- EOT))
# plotted in grid-index space (i from 70ᵒE eastward, j south to north), rows with φ < 64ᵒ (the lat-lon part of the tripolar grid)
J = findall(j -> φ[1, j] < 64, 1:size(φ, 2))
fig = Figure(size = (1900, 1150), fontsize = 13)
top = GridLayout(fig[1, 1:length(cols)])
ax = Axis(top[1, 1]; width = 640, height = 300, title = "EOT20 M2 amplitude (truth)", titlesize = 12, xticksvisible = false, yticksvisible = false,
          xticklabelsvisible = false, yticklabelsvisible = false)
truth = heatmap!(ax, 1:size(EOT, 1), J, ifelse.((depth .> 0)[:, J], abs.(EOT[:, J]), NaN); colormap = :viridis, colorrange = (0, 1.5), nan_color = :gray80)
Colorbar(top[1, 2], truth; height = 300, label = "EOT20 M2 amplitude (m)")
for (r, rowname) in enumerate(rows), (c, colname) in enumerate(cols)
    k = findfirst(i -> runs[i, 1] == rowname && runs[i, 2] == colname, 1:size(runs, 1))
    isnothing(k) && continue
    name = runs[k, 3]
    s = get(scores, name, fill(NaN, 20))
    ax = Axis(fig[r + 1, c]; title = @sprintf("%s, %s\nEOT20 deep %.2f, TPXO9 deep %.2f, gauges %.1f / %.1f cm", rowname, colname,
                                          s[col("EOT_deep")], s[col("TPXO_deep")], s[col("gauge_all")], s[col("gauge_200")]),
              titlesize = 12, xticksvisible = false, yticksvisible = false, xticklabelsvisible = false, yticklabelsvisible = false)
    A = isfile(joinpath(output, name * "_harmonics.jld2")) ? error_magnitude(name) : fill(NaN, size(depth))
    wet = (r == 1 ? mom6depth : depth) .> 0
    hm = heatmap!(ax, 1:size(A, 1), J, ifelse.(wet[:, J], A[:, J], NaN); colormap = :inferno, colorrange = (0, 0.5), highclip = :white, nan_color = :gray80)
    (r, c) == (1, 1) && Colorbar(fig[2:3, length(cols) + 1], hm; label = "|Z_model − Z_EOT20| (m), complex M2 error")
end
Label(fig[0, 1:length(cols)], "EOT20 M2 amplitude and complex M2 error |Z − Z_EOT20| (amplitude and phase), single-layer 0.36ᵒ model (fit days 30–60); errors rms|ΔZ|/√2: deep (>1000 m, |φ|<66ᵒ) vs EOT20 and TPXO9, TICON-3 gauges all / >200 m"; fontsize = 15)
save(joinpath(here, "matrix_amplitude.png"), fig; px_per_unit = 1.5)

fig = Figure(size = (1500, 900), fontsize = 13)
labels = [@sprintf("%s\n%s", runs[i, 1], runs[i, 2]) for i in 1:size(runs, 1)]
names = runs[:, 3]
metrics = [("EOT_deep", "EOT20 deep"), ("TPXO_deep", "TPXO9 deep"), ("EOT_global", "EOT20 global"), ("gauge_200", "gauges > 200 m")]
ax = Axis(fig[1, 1]; limits = (nothing, (0, 38)), ylabel = "rms |ΔZ|/√2 (cm)", title = "M2 errors", xticks = (1:length(names), labels), xticklabelrotation = 0, xticklabelsize = 10)
colors = Makie.wong_colors()
for (q, (m, l)) in enumerate(metrics)
    vals = [get(scores, n, fill(NaN, 20))[col(m)] for n in names]
    barplot!(ax, (1:length(names)) .+ (q - 2.5) * 0.2, vals; width = 0.19, color = colors[q], label = l)
end
axislegend(ax; position = :lt, nbanks = 2)
ax2 = Axis(fig[2, 1]; ylabel = "wall s per model day (T4)", yscale = log10, title = "Cost on one T4: whole run and SAL overhead (≤ 0.1 shown at 0.1)",
           xticks = (1:length(names), labels), xticklabelsize = 10)
tot = [get(costs, n, fill(NaN, 4))[1] for n in names]
none_of(n) = names[findfirst(i -> runs[i, 1] == runs[findfirst(==(n), names), 1] && runs[i, 2] == "none", 1:length(names))]
sal = [get(costs, n, fill(NaN, 4))[1] - get(costs, none_of(n), fill(NaN, 4))[1] for n in names]   # overhead vs no SAL
barplot!(ax2, (1:length(names)) .- 0.2, tot; fillto = 0.1, width = 0.38, color = colors[1], label = "whole run")
barplot!(ax2, (1:length(names)) .+ 0.2, max.(ifelse.(isfinite.(sal), sal, 0.1), 0.1); fillto = 0.1, width = 0.38, color = colors[2], label = "SAL overhead (minus no-SAL run)")
axislegend(ax2; position = :lt)
save(joinpath(here, "matrix_bars.png"), fig; px_per_unit = 1.5)
