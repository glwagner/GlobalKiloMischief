# M2 from the global tide runs compared with tide gauges, with EOT20 as a baseline.
#
#   julia --project=sal_env analysis/m2_gauges.jl
#
# Gauge sets
#   TICON-3 (egbert_ray/data/TICON3/TICON_3.txt): tab-separated; columns lat, lon, constituent, amplitude (cm),
#     Greenwich phase lag (deg, [0, 360)), σA, σG, % missing, N obs, max gap (d), first date, last date, source,
#     gauge type (Coastal | River | Lake). Only "Coastal" gauges are used.
#   Ray (2013) pelagic bottom-pressure set (egbert_ray/data/Ray2013_BPR): read if a file is there (whitespace
#     columns; BPR_COLUMNS="lat,lon,amp_cm,phase_deg" gives their positions, default 1,2,3,4; lines that do not parse
#     are skipped). Same convention: amplitude in cm, Greenwich phase lag in degrees.
#
# Models: M2 is fitted (with S2, K1, O1 and a mean) to hourly η over days 5–20 exactly as in barotropic_m2.jl, and
# converted to A exp(-iG) with G = V₀ + u − arg Z and A = |Z| / f. The 1/8ᵒ stratified runs use the complex M2 η
# amplitudes already fitted by sal_harmonics.jl (analysis/sal_harmonics.jld2, same fit, same days).
# Each gauge is matched to the nearest wet cell of each model grid; gauges with no wet cell within 1.5 cells are
# dropped. EOT20 is interpolated bilinearly (nearest finite of the four neighbours if one is missing).
#
# Statistics (on the set of gauges valid for every model and EOT20, so rows are comparable):
#   RMS misfit  = sqrt(mean |Z_m − Z_o|² / 2)  (the RSS/√2 convention of Ray 2013 / Stammer et al. 2014: the rms of
#                 the difference of the two time series), with Z = A exp(-iG)
#   ΔA          = mean and rms of A_m − A_o
#   phase offset= arg Σ Z_m conj(Z_o) (amplitude-weighted mean phase difference; positive = model lags).
# Prints a table and saves analysis/m2_gauges.png and analysis/m2_gauges.jld2 (gauge-level values, small).

using NumericalEarth
using Oceananigans
using Oceananigans.Grids: λnode, φnode
using JLD2, CairoMakie, Printf, Statistics, Dates, LinearAlgebra, DelimitedFiles

const NCDatasets = Base.require(Base.PkgId(Base.UUID("85f8d34a-cbdd-5861-8df4-14fed0d494ab"), "NCDatasets"))
const NN = Base.require(Base.PkgId(Base.UUID("b8a86587-4115-5ab1-83bc-aa920d37bbce"), "NearestNeighbors"))

first_day, last_day = 5.0, 20.0
here = @__DIR__
data = joinpath(here, "..", "egbert_ray", "data")
output = joinpath(here, "..", "GlobalKiloMischief", "tides", "output")

#####
##### Gauges
#####

ticon = readdlm(joinpath(data, "TICON3", "TICON_3.txt"), '\t', Any)
rows = [r for r in eachrow(ticon) if r[3] == "M2" && strip(string(r[14])) == "Coastal"]
number(x) = x isa Number ? Float64(x) : something(tryparse(Float64, strip(string(x))), NaN)
ticon_gauges = (φ = [number(r[1]) for r in rows], λ = [mod(number(r[2]), 360) for r in rows],
                Z = [number(r[4]) / 100 * cis(-deg2rad(number(r[5]))) for r in rows])
@info "TICON-3: $(length(rows)) coastal M2 gauges"

function read_bpr()
    files = filter(f -> !startswith(f, "."), readdir(joinpath(data, "Ray2013_BPR"); join = true))
    columns = parse.(Int, split(get(ENV, "BPR_COLUMNS", "1,2,3,4"), ","))
    φ, λ, Z = Float64[], Float64[], ComplexF64[]
    for file in files, line in eachline(file)
        fields = split(strip(line))
        length(fields) ≥ maximum(columns) || continue
        values = tryparse.(Float64, fields[columns])
        any(isnothing, values) && continue
        push!(φ, values[1]); push!(λ, mod(values[2], 360)); push!(Z, values[3] / 100 * cis(-deg2rad(values[4])))
    end
    return (; φ, λ, Z)
end
bpr = read_bpr()
@info "Ray (2013) BPR: $(length(bpr.Z)) gauges" * (isempty(bpr.Z) ? " (directory empty: pelagic set skipped)" : "")

gauge_sets = [("TICON", ticon_gauges)]
isempty(bpr.Z) || push!(gauge_sets, ("BPR", bpr))

#####
##### Grids, nearest wet cell
#####

unit(λ, φ) = [cosd(φ) * cosd(λ), cosd(φ) * sind(λ), sind(φ)]

grids = Dict{Float64, Any}()
function grid_nodes(cells_per_degree)
    get!(grids, cells_per_degree) do
        Nx, Ny = round(Int, 360 * cells_per_degree), round(Int, 170 * cells_per_degree)
        grid = TripolarGrid(CPU(); size = (Nx, Ny, 1), halo = (7, 7, 4), z = (-6000, 0))
        λ = [mod(λnode(i, j, 1, grid, Center(), Center(), Center()), 360) for i in 1:Nx, j in 1:Ny]
        φ = [φnode(i, j, 1, grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
        (; grid, λ, φ, Nx, Ny)
    end
end

"""Linear index of the nearest wet cell to each gauge, or 0 if none within 1.5 cells."""
function nearest_wet(g, wet, cells_per_degree, gauges)
    cells = findall(vec(wet))
    points = reduce(hcat, [unit(g.λ[c], g.φ[c]) for c in cells])
    tree = NN.KDTree(points)
    targets = reduce(hcat, [unit(λ, φ) for (λ, φ) in zip(gauges.λ, gauges.φ)])
    idx, dist = NN.nn(tree, targets)
    limit = 2sind(1.5 / cells_per_degree / 2) # chord length of 1.5 cells (in latitude)
    return [d ≤ limit ? cells[i] : 0 for (i, d) in zip(idx, dist)]
end

#####
##### Harmonic fit and Greenwich phase (as in barotropic_m2.jl)
#####

periods = (M2 = 12.4206012, S2 = 12.0, K1 = 23.9344696, O1 = 25.8193417) # hours
ω = [2π / (3600T) for T in values(periods)]

function m2_fit(path, Nx, Ny)
    jldopen(path) do file
        iterations = keys(file["timeseries/t"])
        times = [file["timeseries/t/$i"] for i in iterations]
        window = findall(t -> first_day * 86400 ≤ t ≤ last_day * 86400, times)
        τ = times[window]
        A = hcat(ones(length(τ)), [cos(ωₖ * t) for t in τ, ωₖ in ω], [sin(ωₖ * t) for t in τ, ωₖ in ω])
        P = pinv(A)
        K = length(ω)
        a = zeros(Nx, Ny); b = zeros(Nx, Ny)
        for (n, i) in enumerate(iterations[window])
            η = Float64.(file["timeseries/η/$i"][8:Nx+7, 8:Ny+7, 1])
            a .+= P[2, n] .* η
            b .+= P[2 + K, n] .* η
        end
        @info @sprintf("%s: fitted %d snapshots, days %.2f–%.2f", basename(path), length(window), τ[1] / 86400, τ[end] / 86400)
        return a .- im .* b
    end
end

julian_centuries(date) = (Dates.datetime2julian(date) - 2451545.0) / 36525
T₀ = julian_centuries(DateTime(2015, 1, 1))
V₀ = mod(2(280.4664567 + 36000.76983T₀) - 2(218.3164477 + 481267.88123421T₀), 360)
N = 125.04452 - 1934.136261julian_centuries(DateTime(2015, 1, 1) + Millisecond(round(Int, (first_day + last_day) / 2 * 86400e3)))
u = -2.14sind(N)
f = 1.0004 - 0.0373cosd(N) + 0.0002cosd(2N)
@info @sprintf("V₀ = %.2fᵒ, u = %.2fᵒ, f = %.4f", V₀, u, f)
greenwich(Z) = abs(Z) / f * cis(angle(Z) - deg2rad(V₀ + u))

#####
##### Runs
#####

barotropic(cpd, name) = joinpath(output, @sprintf("barotropic_tides_%gdeg_dt10min_%s.jld2", 1 / cpd, name))
c36 = 1 / 0.36
runs = Any[]
for (sal, label) in [("none", "none"), ("scalar", "scalar"), ("spherical_harmonic", "SH n=40"),
                     ("convolution_prem", "conv PREM"), ("convolution_paper", "conv paper")]
    push!(runs, ("0.36° 1L $label", c36, barotropic(c36, sal)))
    push!(runs, ("0.36° 1L $label χ=3", c36, barotropic(c36, sal * "_wavedrag3")))
end
push!(runs, ("0.25° 1L SH n=40 χ=3", 4.0, barotropic(4.0, "spherical_harmonic_wavedrag3")))
stratified = [("1/8° 10L none", "internal_tides_CUDA_Float64_0.125deg_Nz10_dt3min_tau30d_catke"),
              ("1/8° 10L scalar", "internal_tides_CUDA_Float64_0.125deg_Nz10_dt3min_tau30d_scalar0.1_catke"),
              ("1/8° 10L SH n=40", "internal_tides_CUDA_Float64_0.125deg_Nz10_dt3min_tau30d_harmonic40_catke")]
filter!(r -> isfile(r[3]), runs)

# Model Greenwich Z at gauges: Dict(run label => Dict(set name => Vector (NaN where no wet cell)))
at_gauges = Dict{String, Dict{String, Vector{ComplexF64}}}()
depth_at = Dict{String, Vector{Float64}}() # 0.36ᵒ model depth at nearest wet cell, for subsets

function sample!(label, cpd, Z)
    g = grid_nodes(cpd)
    wet = abs.(Z) .> 0
    at_gauges[label] = Dict(name => [k == 0 ? complex(NaN) : greenwich(Z[k]) for k in nearest_wet(g, wet, cpd, gauges)]
                            for (name, gauges) in gauge_sets)
end

for (label, cpd, path) in runs
    g = grid_nodes(cpd)
    sample!(label, cpd, m2_fit(path, g.Nx, g.Ny))
end

saved = jldopen(joinpath(here, "sal_harmonics.jld2"))
@assert (saved["first_day"], saved["last_day"]) == (first_day, last_day) "sal_harmonics.jld2 fitted a different window"
amplitudes = saved["amplitudes"]
close(saved)
for (label, run) in stratified
    key = "$run/η"
    haskey(amplitudes, key) || (@warn "no $key in sal_harmonics.jld2"; continue)
    sample!(label, 8.0, amplitudes[key]["M2"])
end
amplitudes = nothing; GC.gc()
model_labels = [[r[1] for r in runs]; [s[1] for s in stratified if haskey(at_gauges, s[1])]]

# Depth from the 0.36ᵒ bathymetry used by the single-layer runs
let g = grid_nodes(c36)
    depth = -interior(regrid_bathymetry(g.grid; minimum_depth = 20, interpolation_passes = 3, major_basins = 1), :, :, 1)
    for (name, gauges) in gauge_sets
        k = nearest_wet(g, depth .> 0, c36, gauges)
        depth_at[name] = [i == 0 ? NaN : depth[i] for i in k]
    end
end

#####
##### EOT20 at gauges
#####

eot = NCDatasets.NCDataset(joinpath(data, "EOT20", "ocean_tides", "M2_ocean_eot20.nc"))
eot_λ = Float64.(eot["lon"][:]); eot_φ = Float64.(eot["lat"][:])
eot_Z = map((A, G) -> ismissing(A) || ismissing(G) ? complex(NaN) : A / 100 * cis(-deg2rad(G)), eot["amplitude"][:, :], eot["phase"][:, :])
close(eot)
Δλ, Δφ = eot_λ[2] - eot_λ[1], eot_φ[2] - eot_φ[1]
function eot_interpolate(λ, φ)
    x = mod(λ - eot_λ[1], 360) / Δλ + 1; y = (φ - eot_φ[1]) / Δφ + 1
    i = clamp(floor(Int, x), 1, length(eot_λ) - 1); j = clamp(floor(Int, y), 1, length(eot_φ) - 1)
    a, b = x - i, y - j
    corners = (eot_Z[i, j], eot_Z[i+1, j], eot_Z[i, j+1], eot_Z[i+1, j+1])
    weights = ((1 - a) * (1 - b), a * (1 - b), (1 - a) * b, a * b)
    all(isfinite, corners) && return sum(weights .* corners)
    finite = findall(isfinite, corners)
    isempty(finite) && return complex(NaN)
    return corners[finite[argmax(collect(weights[finite]))]]
end
for (name, gauges) in gauge_sets
    get!(at_gauges, "EOT20", Dict{String, Vector{ComplexF64}}())[name] = eot_interpolate.(gauges.λ, gauges.φ)
end
all_labels = ["EOT20"; model_labels]

#####
##### Statistics
#####

function statistics(Zm, Zo)
    D = Zm .- Zo
    ΔA = abs.(Zm) .- abs.(Zo)
    (n = length(D), rms = 100sqrt(mean(abs2, D) / 2), mean_ΔA = 100mean(ΔA), rms_ΔA = 100sqrt(mean(abs2, ΔA)),
     phase = rad2deg(angle(sum(Zm .* conj.(Zo)))), rms_signal = 100sqrt(mean(abs2, Zo) / 2))
end

subsets = Dict{String, Any}()
results = Dict{Tuple{String, String}, Any}()
for (name, gauges) in gauge_sets
    valid = trues(length(gauges.Z))
    for label in all_labels
        valid .&= isfinite.(at_gauges[label][name])
    end
    valid .&= isfinite.(gauges.Z) .& (abs.(gauges.Z) .> 0)
    d = depth_at[name]
    sets = name == "TICON" ?
        [("TICON all", valid), ("TICON open (depth₀.₃₆ > 200 m)", valid .& (d .> 200)),
         ("TICON deep (depth₀.₃₆ > 1000 m)", valid .& (d .> 1000))] :
        [("BPR (Ray 2013)", valid)]
    for (setname, mask) in sets
        subsets[setname] = (name, mask)
        for label in all_labels
            results[(label, setname)] = statistics(at_gauges[label][name][mask], gauges.Z[mask])
        end
    end
end

setnames = sort(collect(keys(subsets)))
println()
println("M2 vs gauges. RMS = sqrt(mean|Z_m − Z_o|²/2) [cm]; ΔA = model − gauge amplitude [cm]; phase = arg Σ Z_m Z_o* [deg, + = model lags]")
for setname in setnames
    first_result = results[(first(all_labels), setname)]
    @printf("\n%s: %d gauges (common to all models), rms signal %.1f cm\n", setname, first_result.n, first_result.rms_signal)
    @printf("| %-24s | %8s | %8s | %8s | %7s |\n", "run", "RMS", "mean ΔA", "rms ΔA", "phase")
    @printf("|%s|%s|%s|%s|%s|\n", "-"^26, "-"^10, "-"^10, "-"^10, "-"^9)
    for label in all_labels
        r = results[(label, setname)]
        @printf("| %-24s | %8.2f | %8.2f | %8.2f | %7.1f |\n", label, r.rms, r.mean_ΔA, r.rms_ΔA, r.phase)
    end
end

jldsave(joinpath(here, "m2_gauges.jld2"); all_labels, at_gauges, depth_at,
        gauges = Dict(name => g for (name, g) in gauge_sets), results, subsets = Dict(k => v[2] for (k, v) in subsets))

#####
##### Figure: maps of |ΔZ|/√2 at TICON gauges and amplitude scatter
#####

name, gauges = gauge_sets[1]
mask = subsets["TICON all"][2]
show = filter(l -> haskey(at_gauges, l), ["EOT20", "0.36° 1L scalar χ=3", "0.36° 1L SH n=40 χ=3", "0.25° 1L SH n=40 χ=3",
                                          "1/8° 10L scalar", "1/8° 10L SH n=40"])
fig = Figure(size = (2100, 900 + 0 * length(show)))
Label(fig[0, 1:length(show)], "M2 at $(count(mask)) TICON-3 coastal gauges: misfit |Z_model − Z_gauge|/√2 (top) and amplitude (bottom)"; fontsize = 22)
for (col, label) in enumerate(show)
    Zm = at_gauges[label][name][mask]; Zo = gauges.Z[mask]
    misfit = 100abs.(Zm .- Zo) ./ √2
    r = results[(label, "TICON all")]
    ax = Axis(fig[1, col]; title = @sprintf("%s\nRMS %.1f cm", label, r.rms), aspect = DataAspect(), limits = (0, 360, -80, 90))
    order = sortperm(misfit)
    sc = scatter!(ax, gauges.λ[mask][order], gauges.φ[mask][order]; color = misfit[order], colormap = :inferno,
                  colorrange = (0, 50), markersize = 5)
    col == length(show) && Colorbar(fig[1, col + 1], sc; label = "|ΔZ|/√2 (cm)")
    ax = Axis(fig[2, col]; xlabel = "gauge amplitude (cm)", ylabel = col == 1 ? "model amplitude (cm)" : "",
              aspect = 1, limits = (0, 400, 0, 400))
    scatter!(ax, 100abs.(Zo), 100abs.(Zm); markersize = 4, color = (:steelblue, 0.5))
    lines!(ax, [0, 400], [0, 400]; color = :black)
end
save(joinpath(here, "m2_gauges.png"), fig)
@info "Saved analysis/m2_gauges.png and analysis/m2_gauges.jld2"
