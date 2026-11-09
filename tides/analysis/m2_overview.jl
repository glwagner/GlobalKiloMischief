# M2 overview for the single-layer runs with true depths (60 days, online harmonic fit over days 30–60):
# top row EOT20 amplitude (left) and model vs tide-gauge amplitude at TICON-3 gauges (right),
# bottom row the complex M2 error |Z_model − Z_EOT20| (amplitude and phase) with no, scalar, spherical harmonic and
# convolution SAL.
# Errors are in the units of Chen et al. (2026), rms |Z_model − Z_obs| / √2.
#   julia --project=sal_env analysis/m2_overview.jl
# RUN_SUFFIX (default _wavedrag2_fix60d_harmonics) selects the runs.

using NumericalEarth
using Oceananigans
using Oceananigans.Grids: λnode, φnode
using Oceananigans.Operators: Azᶜᶜᶜ
using JLD2, CairoMakie, Printf, Statistics, Dates, DelimitedFiles

const NCDatasets = Base.require(Base.PkgId(Base.UUID("85f8d34a-cbdd-5861-8df4-14fed0d494ab"), "NCDatasets"))
const NN = Base.require(Base.PkgId(Base.UUID("b8a86587-4115-5ab1-83bc-aa920d37bbce"), "NearestNeighbors"))

here = @__DIR__
data = joinpath(here, "..", "egbert_ray", "data")
output = joinpath(here, "..", "GlobalKiloMischief", "tides", "output")
suffix = get(ENV, "RUN_SUFFIX", "_wavedrag2_fix60d_harmonics")
runs = [("no SAL", "none"), ("scalar SAL (β = 0.1)", "scalar"), ("spherical harmonic SAL (n = 40)", "spherical_harmonic"),
        ("convolution SAL (PREM)", "convolution_prem")]
run_file(sal) = joinpath(output, "barotropic_tides_0.36deg_dt10min_$(sal)$(suffix).jld2")
# OVERVIEW_RUNS = four run names (without _harmonics.jld2), comma separated, in the order of `runs`, overrides RUN_SUFFIX;
# OVERVIEW_CONFIG names the configuration in the title.
if haskey(ENV, "OVERVIEW_RUNS")
    let names = split(ENV["OVERVIEW_RUNS"], ","), labels = [r[1] for r in runs]
        global runs = [(l, String(n)) for (l, n) in zip(labels, names)]
        global run_file(name) = joinpath(output, name * "_harmonics.jld2")
    end
end
config = get(ENV, "OVERVIEW_CONFIG", "wave drag χ = 2")
filter!(r -> isfile(run_file(r[2])) || (@warn "missing $(run_file(r[2]))"; false), runs)

cells_per_degree = 1 / 0.36
Nx, Ny = round(Int, 360 * cells_per_degree), round(Int, 170 * cells_per_degree)
grid = TripolarGrid(CPU(); size = (Nx, Ny, 1), halo = (7, 7, 4), z = (-6000, 0))
depth = -interior(regrid_bathymetry(grid; minimum_depth = 20, interpolation_passes = 3, major_basins = 1), :, :, 1)
wet = depth .> 0
λ = [mod(λnode(i, j, 1, grid, Center(), Center(), Center()), 360) for i in 1:Nx, j in 1:Ny]
φ = [φnode(i, j, 1, grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
area = [Azᶜᶜᶜ(i, j, 1, grid) for i in 1:Nx, j in 1:Ny]

##### EOT20 at cell centres

eot = NCDatasets.NCDataset(joinpath(data, "EOT20", "ocean_tides", "M2_ocean_eot20.nc"))
eot_λ = Float64.(eot["lon"][:]); eot_φ = Float64.(eot["lat"][:])
eot_Z = map((A, G) -> ismissing(A) || ismissing(G) ? complex(NaN) : A / 100 * cis(-deg2rad(G)), eot["amplitude"][:, :], eot["phase"][:, :])
close(eot)
Δλ, Δφ = eot_λ[2] - eot_λ[1], eot_φ[2] - eot_φ[1]
function interpolate(λ, φ)
    x = mod(λ - eot_λ[1], 360) / Δλ + 1; y = (φ - eot_φ[1]) / Δφ + 1
    i = clamp(floor(Int, x), 1, length(eot_λ) - 1); j = clamp(floor(Int, y), 1, length(eot_φ) - 1)
    a, b = x - i, y - j
    return (1 - a) * (1 - b) * eot_Z[i, j] + a * (1 - b) * eot_Z[i+1, j] + (1 - a) * b * eot_Z[i, j+1] + a * b * eot_Z[i+1, j+1]
end
observed = interpolate.(λ, φ)
compared = wet .& isfinite.(observed)
deep = compared .& (depth .> 1000) .& (abs.(φ) .< 66)
rms(x, mask) = sqrt(sum(abs2.(x[mask]) .* area[mask]) / sum(area[mask]))

##### Model M2 in Greenwich phase (as in diag/m2_table.jl)

julian_centuries(date) = (Dates.datetime2julian(date) - 2451545.0) / 36525
T₀ = julian_centuries(DateTime(2015, 1, 1))
V₀ = mod(2(280.4664567 + 36000.76983T₀) - 2(218.3164477 + 481267.88123421T₀), 360)
function greenwich(Z, window_mid)
    N = 125.04452 - 1934.136261julian_centuries(DateTime(2015, 1, 1) + Millisecond(round(Int, window_mid * 86400e3)))
    u = -2.14sind(N); f = 1.0004 - 0.0373cosd(N) + 0.0002cosd(2N)
    return @. abs(Z) / f * cis(angle(Z) - deg2rad(V₀ + u))
end
models = map(runs) do (label, sal)
    h = load(run_file(sal))
    @info @sprintf("%s: %d samples, days %.1f–%.1f", label, h["samples"], h["first_day"], h["last_day"])
    greenwich(h["Z"]["M2"], (h["first_day"] + h["last_day"]) / 2) .* wet
end

##### TICON-3 coastal gauges, matched to the nearest wet cell within 1.5 cells

ticon = readdlm(joinpath(data, "TICON3", "TICON_3.txt"), '\t', Any)
rows = [r for r in eachrow(ticon) if r[3] == "M2" && strip(string(r[14])) == "Coastal"]
number(x) = x isa Number ? Float64(x) : something(tryparse(Float64, strip(string(x))), NaN)
gφ = [number(r[1]) for r in rows]; gλ = [mod(number(r[2]), 360) for r in rows]
gZ = [number(r[4]) / 100 * cis(-deg2rad(number(r[5]))) for r in rows]
unit(λ, φ) = [cosd(φ) * cosd(λ), cosd(φ) * sind(λ), sind(φ)]
cells = findall(vec(wet))
tree = NN.KDTree(reduce(hcat, [unit(λ[c], φ[c]) for c in cells]))
idx, dist = NN.nn(tree, reduce(hcat, [unit(a, b) for (a, b) in zip(gλ, gφ)]))
nearest = [d ≤ 2sind(1.5 * 0.36 / 2) ? cells[i] : 0 for (i, d) in zip(idx, dist)]
valid = (nearest .> 0) .& isfinite.(gZ) .& (abs.(gZ) .> 0)
valid .&= [n > 0 && all(m -> isfinite(m[n]), models) && isfinite(observed[n]) for n in nearest]
open_coast = valid .& [n > 0 && depth[n] > 200 for n in nearest]
gauge_rms(field, mask) = 100sqrt(mean(abs2, field[nearest[mask]] .- gZ[mask]) / 2)
@info "TICON-3: $(count(valid)) gauges matched, $(count(open_coast)) where the model depth exceeds 200 m"

##### Figure

masked(x, mask) = ifelse.(mask, x, NaN)
deep_error(m) = 100rms(m .- observed, deep) / √2
colorrange = (0, 1.2)
# Top row: EOT20 map with its colorbar, then the gauge scatter at the map's height.
# Bottom block: the four error maps in a 2 × 2 grid with one shared colorbar. The configuration goes in the caption.
fig = Figure(size = (2000, 1650), fontsize = 20)
top = fig[1, 1] = GridLayout()
bottom = fig[2, 1] = GridLayout()

ax = Axis(top[1, 1]; title = @sprintf("EOT20 M2 amplitude (deep rms %.1f cm) and TICON-3 gauges", 100rms(abs.(observed), deep)),
          aspect = DataAspect())
hidedecorations!(ax)
hm = heatmap!(ax, masked(abs.(observed), compared); colormap = :viridis, colorrange)
scatter!(ax, [i for i in 1:Nx, j in 1:Ny][nearest[valid]], [j for i in 1:Nx, j in 1:Ny][nearest[valid]];
         markersize = 6, color = (:white, 0.8), strokecolor = :black, strokewidth = 0.5)
Colorbar(top[1, 2], hm; label = "amplitude (m)", height = Relative(0.7))

ax = Axis(top[1, 3]; title = "Model vs gauge amplitude, $(count(open_coast)) gauges deeper than 200 m",
          xlabel = "gauge amplitude (m)", ylabel = "model amplitude (m)", aspect = 1)
lim = 1.05maximum(abs, gZ[open_coast])
lines!(ax, [0, lim], [0, lim]; color = :black, linestyle = :dash)
colors = Makie.wong_colors()
scatter!(ax, abs.(gZ[open_coast]), abs.(observed[nearest[open_coast]]); color = (:gray, 0.6), markersize = 7,
         label = @sprintf("EOT20: %.1f cm", gauge_rms(observed, open_coast)))
for (k, ((label, _), m)) in enumerate(zip(runs, models))
    scatter!(ax, abs.(gZ[open_coast]), abs.(m[nearest[open_coast]]); color = (colors[k], 0.6), markersize = 7,
             label = @sprintf("%s: %.1f cm", label, gauge_rms(m, open_coast)))
end
limits!(ax, 0, lim, 0, lim)
axislegend(ax, "rms |Z − Z_gauge|/√2"; position = :rb, framevisible = false, labelsize = 16, titlesize = 16)
colsize!(top, 1, Auto(2.2)); colsize!(top, 3, Auto(1))

for (k, ((label, _), m)) in enumerate(zip(runs, models))
    local ax = Axis(bottom[fld1(k, 2), mod1(k, 2)]; aspect = DataAspect(),
                    title = @sprintf("%s: deep %.2f cm, all gauges %.1f cm", label, deep_error(m), gauge_rms(m, valid)))
    hidedecorations!(ax)
    global err = heatmap!(ax, masked(abs.(m .- observed), compared); colormap = :inferno, colorrange = (0, 0.5), highclip = :white)
end
Colorbar(bottom[1:2, 3], err; label = "|Z_model − Z_EOT20| (m)", height = Relative(0.8))
colgap!(top, 15); colgap!(bottom, 15); rowgap!(bottom, 15); rowgap!(fig.layout, 25)
resize_to_layout!(fig)
save(joinpath(here, get(ENV, "OVERVIEW_PNG", "m2_overview.png")), fig)
@info "Saved analysis/m2_overview.png"
for ((label, _), m) in zip(runs, models)
    println(@sprintf("%-34s deep vs EOT20 %5.2f cm | gauges all %6.2f cm, depth > 200 m %6.2f cm", label, deep_error(m),
                     gauge_rms(m, valid), gauge_rms(m, open_coast)))
end
