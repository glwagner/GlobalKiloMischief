# M2 of finished *_harmonics runs against EOT20, GOT5.5, HAMTIDE11a (GOT5.6 holds only third-degree constituents; TPXO
# needs registration). Errors in paper units rms|ΔZ|/√2 (cm): plain area mean over the mask and literal ÷4πR² (eq. 22).
#   CELLS_PER_DEGREE=2.777.. julia --project=sal_env analysis/diag/m2_atlases.jl <run>...
# Appends rows to stdout as TSV. All scores on the mask common to the three atlases, and each atlas's own.
using NumericalEarth, Oceananigans
using Oceananigans.Grids: λnode, φnode
using Oceananigans.Operators: Azᶜᶜᶜ
using JLD2, Printf, Dates
const NCDatasets = Base.require(Base.PkgId(Base.UUID("85f8d34a-cbdd-5861-8df4-14fed0d494ab"), "NCDatasets"))
output = joinpath(@__DIR__, "..", "..", "GlobalKiloMischief", "tides", "output")
data = joinpath(@__DIR__, "..", "..", "egbert_ray", "data")
cpd = parse(Float64, get(ENV, "CELLS_PER_DEGREE", "2.7777777777777777"))
Nx, Ny = round(Int, 360cpd), round(Int, 170cpd)
gridfile = joinpath(@__DIR__, @sprintf("grid_%d.jld2", Nx))
if isfile(gridfile)
    depth, λ, φ, area = load(gridfile, "depth", "λ", "φ", "area")
else
    grid = TripolarGrid(CPU(); size = (Nx, Ny, 1), halo = (7, 7, 4), z = (-6000, 0))
    depth = -interior(regrid_bathymetry(grid; minimum_depth = 20, interpolation_passes = 3, major_basins = 1), :, :, 1)
    λ = [mod(λnode(i, j, 1, grid, Center(), Center(), Center()), 360) for i in 1:Nx, j in 1:Ny]
    φ = [φnode(i, j, 1, grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
    area = [Azᶜᶜᶜ(i, j, 1, grid) for i in 1:Nx, j in 1:Ny]
    jldsave(gridfile; depth, λ, φ, area)
end
wet = depth .> 0
valid(x) = !ismissing(x) && isfinite(x) && abs(x) < 1e6
function atlas(path, lon, lat, re, im; polar = false)
    ds = NCDatasets.NCDataset(path)
    lo = Float64.(ds[lon][:]); la = Float64.(ds[lat][:])
    return lo, la, ds
end
function load_atlas(name)
    if name == "EOT20"
        ds = NCDatasets.NCDataset(joinpath(data, "EOT20/ocean_tides/M2_ocean_eot20.nc"))
        lo, la = Float64.(ds["lon"][:]), Float64.(ds["lat"][:]); A, G = ds["amplitude"][:, :], ds["phase"][:, :]
    elseif name == "GOT5.5"
        ds = NCDatasets.NCDataset(joinpath(data, "GOT5.5/GOT5.5/ocean_tides/m2.nc"))
        lo, la = Float64.(ds["longitude"][:]), Float64.(ds["latitude"][:]); A, G = ds["amplitude"][:, :], ds["phase"][:, :]
    else
        ds = NCDatasets.NCDataset(joinpath(data, "HAMTIDE11a/m2.hamtide11a.nc"))
        lo, la = Float64.(ds["LON"][:]), Float64.(ds["LAT"][:]); A, G = ds["AMPL"][:, :], ds["PHAS"][:, :]
    end
    Z = map((a, g) -> valid(a) && valid(g) ? a / 100 * cis(-deg2rad(g)) : complex(NaN), A, G)
    close(ds)
    Δλ, Δφ = lo[2] - lo[1], la[2] - la[1]
    function at(λ, φ)
        x = (λ - lo[1]) / Δλ + 1; y = (φ - la[1]) / Δφ + 1
        i = clamp(floor(Int, x), 1, length(lo) - 1); j = clamp(floor(Int, y), 1, length(la) - 1)
        a, b = x - i, y - j
        return (1 - a) * (1 - b) * Z[i, j] + a * (1 - b) * Z[i+1, j] + (1 - a) * b * Z[i, j+1] + a * b * Z[i+1, j+1]
    end
    return at.(λ, φ)
end
# TPXO9 M2 as regridded by Chen et al. onto MOM6 p04arc9: nearest cell (within 0.5ᵒ) on a 0.1ᵒ lookup
function load_tpxo()
    d = joinpath(@__DIR__, "..", "..", "mom6_repro", "extract", "run_inputs", "p04arc9")
    hg = NCDatasets.NCDataset(joinpath(d, "hgrid_p04arc9.nc")); lo = mod.(Float64.(hg["x"][2:2:end, 2:2:end]), 360); la = Float64.(hg["y"][2:2:end, 2:2:end]); close(hg)
    ds = NCDatasets.NCDataset(joinpath(d, "tpxo9_m2_p04arc9.nc")); A = Float64.(ds["ha"][:, :]); G = Float64.(ds["hp"][:, :]); close(ds)
    @info "TPXO9 ha extrema" extrema(filter(isfinite, A)) extrema(filter(isfinite, G))
    Z = map((a, p) -> isfinite(a) && a > 0 && abs(a) < 1e3 ? a * cis(-deg2rad(p)) : complex(NaN), A, G)
    table = Dict{Tuple{Int, Int}, Vector{Int}}()
    for q in eachindex(lo); push!(get!(table, (floor(Int, lo[q]), floor(Int, la[q])), Int[]), q); end
    function at(x, y)
        best, bd = complex(NaN), 0.5^2
        for a in -1:1, b in -1:1, q in get(table, (mod(floor(Int, x) + a, 360), floor(Int, y) + b), Int[])
            dλ = mod(lo[q] - x + 180, 360) - 180; d = (dλ * cosd(y))^2 + (la[q] - y)^2
            if d < bd; best = Z[q]; bd = d; end
        end
        return best
    end
    return at.(λ, φ)
end
atlases = Dict(n => load_atlas(n) for n in ("EOT20", "GOT5.5", "HAMTIDE"))
get(ENV, "TPXO", "true") == "true" && (atlases["TPXO9"] = load_tpxo())
jc(d) = (Dates.datetime2julian(d) - 2451545.0) / 36525
T₀ = jc(DateTime(2015, 1, 1))
V₀ = mod(2(280.4664567 + 36000.76983T₀) - 2(218.3164477 + 481267.88123421T₀), 360)
function greenwich(Z, mid)
    N = 125.04452 - 1934.136261jc(DateTime(2015, 1, 1) + Millisecond(round(Int, mid * 86400e3)))
    u = -2.14sind(N); f = 1.0004 - 0.0373cosd(N) + 0.0002cosd(2N)
    return @. abs(Z) / f * cis(angle(Z) - deg2rad(V₀ + u))
end
common = wet .& reduce(.&, [isfinite.(a) for a in values(atlases)])
sphere = 4π * 6371e3^2
rms(x, m) = sqrt(sum(abs2.(x[m]) .* area[m]) / sum(area[m]))
f(m) = sqrt(sum(area[m]) / sphere)
println("# area fraction common global ", round(f(common)^2, digits = 3), " deep ", round(f(common .& (depth .> 1000) .& (abs.(φ) .< 66))^2, digits = 3))
println("# atlas mutual (common mask, deep / global, plain area mean, cm)")
ns = collect(keys(atlases)); deep = common .& (depth .> 1000) .& (abs.(φ) .< 66)
for a in ns, b in ns
    a < b && @printf("# %s - %s: %.2f / %.2f\n", a, b, 100rms(atlases[a] .- atlases[b], deep) / √2, 100rms(atlases[a] .- atlases[b], common) / √2)
end
println("run\tatlas\tdeep\tdeep_amp\tdeep_phase\tglobal\tglobal_amp\tglobal_phase\tdeep/4piR2\tglobal/4piR2\tphase_offset_deg\trmsamp_deep")
for name in ARGS
    h = load(joinpath(output, name * ".jld2"))
    m = greenwich(h["Z"]["M2"], (h["first_day"] + h["last_day"]) / 2) .* wet
    short = replace(name, "barotropic_tides_" => "", "dt10min_" => "", "_fix60d_harmonics" => "")
    for n in keys(atlases)
        o = atlases[n]
        t, a = 100rms(m .- o, deep) / √2, 100rms(abs.(m) .- abs.(o), deep) / √2
        tg, ag = 100rms(m .- o, common) / √2, 100rms(abs.(m) .- abs.(o), common) / √2
        off = rad2deg(angle(sum(m[deep] .* conj.(o[deep]) .* area[deep])))
        @printf("%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.1f\t%.2f\n", short, n, t, a, sqrt(max(t^2 - a^2, 0)), tg, ag, sqrt(max(tg^2 - ag^2, 0)),
                t * f(deep), tg * f(common), off, 100rms(abs.(m), deep))
    end
    flush(stdout)
end
