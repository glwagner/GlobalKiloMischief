# Scores M2 of *_harmonics runs vs EOT20 and TPXO9 (deep: depth > 1000 m, |φ| < 66ᵒ; global: all compared cells) and TICON-3
# coastal gauges (all matched, and model depth > 200 m), in rms|ΔZ|/√2 (cm, plain area mean), with the amplitude part
# rms(|m| - |o|)/√2, the remainder sqrt(t² - a²) ("phase"), and the area-weighted mean phase offset arg Σ m o* A (ᵒ).
# Mask: wet in both our regrid_bathymetry depth (grid_1000.jld2) and the MOM6 depth regridded onto our grid, so that runs on
# either bathymetry are scored on the same cells. Depth criteria use our bathymetry.
#   julia --project=sal_env analysis/matrix/score_matrix.jl <run names (without _harmonics.jld2)>...
using JLD2, Printf, Dates, DelimitedFiles, Statistics
const NCDatasets = Base.require(Base.PkgId(Base.UUID("85f8d34a-cbdd-5861-8df4-14fed0d494ab"), "NCDatasets"))
const NN = Base.require(Base.PkgId(Base.UUID("b8a86587-4115-5ab1-83bc-aa920d37bbce"), "NearestNeighbors"))
here = @__DIR__
output = joinpath(here, "..", "..", "GlobalKiloMischief", "tides", "output")
data = joinpath(here, "..", "..", "egbert_ray", "data")
depth, λ, φ, area = load(joinpath(here, "..", "diag", "grid_1000.jld2"), "depth", "λ", "φ", "area")
Nx, Ny = size(depth)
mom6depth = load(joinpath(here, "..", "..", "GlobalKiloMischief", "tides", "mom6_p04arc9_on_1000.jld2"), "depth")
wet = (depth .> 0) .& (mom6depth .> 0)

eot = NCDatasets.NCDataset(joinpath(data, "EOT20/ocean_tides/M2_ocean_eot20.nc"))
lo, la = Float64.(eot["lon"][:]), Float64.(eot["lat"][:])
EZ = map((a, g) -> ismissing(a) || ismissing(g) ? complex(NaN) : a / 100 * cis(-deg2rad(g)), eot["amplitude"][:, :], eot["phase"][:, :]); close(eot)
function bilinear(x0, y0)
    x = mod(x0 - lo[1], 360) / (lo[2] - lo[1]) + 1; y = (y0 - la[1]) / (la[2] - la[1]) + 1
    i = clamp(floor(Int, x), 1, length(lo) - 1); j = clamp(floor(Int, y), 1, length(la) - 1); a, b = x - i, y - j
    (1 - a) * (1 - b) * EZ[i, j] + a * (1 - b) * EZ[i+1, j] + (1 - a) * b * EZ[i, j+1] + a * b * EZ[i+1, j+1]
end
EOT = bilinear.(λ, φ)
function load_tpxo()
    d = joinpath(here, "..", "..", "mom6_repro", "extract", "run_inputs", "p04arc9")
    hg = NCDatasets.NCDataset(joinpath(d, "hgrid_p04arc9.nc")); LO = mod.(Float64.(hg["x"][2:2:end, 2:2:end]), 360); LA = Float64.(hg["y"][2:2:end, 2:2:end]); close(hg)
    ds = NCDatasets.NCDataset(joinpath(d, "tpxo9_m2_p04arc9.nc")); A = Float64.(ds["ha"][:, :]); G = Float64.(ds["hp"][:, :]); close(ds)
    Z = map((a, p) -> isfinite(a) && a > 0 && abs(a) < 1e3 ? a * cis(-deg2rad(p)) : complex(NaN), A, G)
    table = Dict{Tuple{Int, Int}, Vector{Int}}()
    for q in eachindex(LO); push!(get!(table, (floor(Int, LO[q]), floor(Int, LA[q])), Int[]), q); end
    function at(x, y)
        best, bd = complex(NaN), 0.5^2
        for a in -1:1, b in -1:1, q in get(table, (mod(floor(Int, x) + a, 360), floor(Int, y) + b), Int[])
            dλ = mod(LO[q] - x + 180, 360) - 180; dd = (dλ * cosd(y))^2 + (LA[q] - y)^2
            dd < bd && (best = Z[q]; bd = dd)
        end
        best
    end
    at.(λ, φ)
end
TPXO = load_tpxo()
atlases = (EOT20 = EOT, TPXO9 = TPXO)

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
gvalid = (nearest .> 0) .& isfinite.(gZ) .& (abs.(gZ) .> 0)
gdeep = gvalid .& [n > 0 && depth[n] > 200 for n in nearest]

jc(d) = (Dates.datetime2julian(d) - 2451545.0) / 36525
T₀ = jc(DateTime(2015, 1, 1)); V₀ = mod(2(280.4664567 + 36000.76983T₀) - 2(218.3164477 + 481267.88123421T₀), 360)
function greenwich(Z, mid)
    N = 125.04452 - 1934.136261jc(DateTime(2015, 1, 1) + Millisecond(round(Int, mid * 86400e3)))
    u = -2.14sind(N); f = 1.0004 - 0.0373cosd(N) + 0.0002cosd(2N)
    @. abs(Z) / f * cis(angle(Z) - deg2rad(V₀ + u))
end
rms(x, m) = sqrt(sum(abs2.(x[m]) .* area[m]) / sum(area[m]))
masks = Dict(n => wet .& isfinite.(o) for (n, o) in pairs(atlases))
println("# cells: EOT20 global $(count(masks[:EOT20])) deep $(count(masks[:EOT20] .& (depth .> 1000) .& (abs.(φ) .< 66))), TPXO9 global $(count(masks[:TPXO9])); gauges $(count(gvalid)) (>200 m: $(count(gdeep)))")
println("run\tEOT_deep\tEOT_deep_amp\tEOT_deep_phase\tEOT_offset\tEOT_global\tTPXO_deep\tTPXO_deep_amp\tTPXO_deep_phase\tTPXO_offset\tTPXO_global\tgauge_all\tgauge_200\tmax_amp_m")
for name in ARGS
    file = joinpath(output, name * "_harmonics.jld2")
    isfile(file) || (println("$name\tMISSING"); continue)
    h = load(file)
    m = greenwich(h["Z"]["M2"], (h["first_day"] + h["last_day"]) / 2)
    vals = Float64[]
    for (n, o) in pairs(atlases)
        g = masks[n] .& isfinite.(m); d = g .& (depth .> 1000) .& (abs.(φ) .< 66)
        t = 100rms(m .- o, d) / √2; a = 100rms(abs.(m) .- abs.(o), d) / √2
        push!(vals, t, a, sqrt(max(t^2 - a^2, 0)), rad2deg(angle(sum(m[d] .* conj.(o[d]) .* area[d]))), 100rms(m .- o, g) / √2)
    end
    gm(mask) = (k = mask .& [isfinite(m[n]) for n in max.(nearest, 1)]; 100sqrt(mean(abs2, m[nearest[k]] .- gZ[k]) / 2))
    push!(vals, gm(gvalid), gm(gdeep), maximum(abs, m[wet .& isfinite.(m)]))
    println(name, "\t", join([@sprintf("%.2f", v) for v in vals], "\t")); flush(stdout)
end
