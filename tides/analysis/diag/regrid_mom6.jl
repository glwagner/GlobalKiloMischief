# Regrid Chen et al.'s MOM6 p04arc9 depth (depth_p04arc9_aw_moreocean_ocnpct0p5.nc) and JSL_RH (JSLH_p04arc9_aw.nc) onto
# our 0.36ᵒ TripolarGrid by Voronoi binning: every MOM6 h-point goes to our nearest cell centre (both grids are ~0.36ᵒ, so
# each of our cells collects ~1 point, approximately area-conservative). A cell is wet if at least half of its points are
# wet in MOM6 (depth > MASKING_DEPTH = -5 m); its depth and JSL_RH are the means over its wet points (depth is clipped at
# run time to MINIMUM_DEPTH). Cells with no point take the nearest MOM6 point within 0.5ᵒ, else land.
# Writes GlobalKiloMischief/tides/mom6_p04arc9_on_<Nx>.jld2 (depth > 0 wet, 0 land; rh in m s⁻¹).
using JLD2, Printf
const NCDatasets = Base.require(Base.PkgId(Base.UUID("85f8d34a-cbdd-5861-8df4-14fed0d494ab"), "NCDatasets"))
d = joinpath(@__DIR__, "..", "..", "mom6_repro", "extract", "run_inputs", "p04arc9")
λ, φ, ourdepth = load(joinpath(@__DIR__, "grid_1000.jld2"), "λ", "φ", "depth")
Nx, Ny = size(λ)
hg = NCDatasets.NCDataset(joinpath(d, "hgrid_p04arc9.nc")); lo = mod.(Float64.(hg["x"][2:2:end, 2:2:end]), 360); la = Float64.(hg["y"][2:2:end, 2:2:end]); close(hg)
ds = NCDatasets.NCDataset(joinpath(d, "depth_p04arc9_aw_moreocean_ocnpct0p5.nc")); D = Float64.(ds["depth"][:, :]); close(ds)
ds = NCDatasets.NCDataset(joinpath(d, "JSLH_p04arc9_aw.nc")); R = Float64.(coalesce.(ds["JSL_RH"][:, :], 0.0)); close(ds)
@assert size(D) == size(lo) "$(size(D)) vs $(size(lo))"
mwet = D .> -5
table = Dict{Tuple{Int, Int}, Vector{Int}}()
for q in eachindex(λ); push!(get!(table, (floor(Int, λ[q]), floor(Int, φ[q])), Int[]), q); end
function nearest(x, y, table, LO, LA, r2)
    best, bd = 0, r2
    for a in -1:1, b in -1:1, q in get(table, (mod(floor(Int, x) + a, 360), floor(Int, y) + b), Int[])
        dλ = mod(LO[q] - x + 180, 360) - 180; dd = (dλ * cosd(y))^2 + (LA[q] - y)^2
        dd < bd && (best = q; bd = dd)
    end
    return best
end
nwet = zeros(Nx, Ny); ntot = zeros(Nx, Ny); sD = zeros(Nx, Ny); sR = zeros(Nx, Ny)
for p in eachindex(lo)
    q = nearest(lo[p], la[p], table, λ, φ, 1.0); q == 0 && continue
    ntot[q] += 1
    mwet[p] && (nwet[q] += 1; sD[q] += D[p]; sR[q] += R[p])
end
mtable = Dict{Tuple{Int, Int}, Vector{Int}}()
for p in eachindex(lo); push!(get!(mtable, (floor(Int, lo[p]), floor(Int, la[p])), Int[]), p); end
depth = zeros(Nx, Ny); rh = zeros(Nx, Ny); filled = 0
for q in eachindex(λ)
    if ntot[q] > 0
        nwet[q] ≥ ntot[q] / 2 && nwet[q] > 0 && (depth[q] = sD[q] / nwet[q]; rh[q] = sR[q] / nwet[q])
    else
        p = nearest(λ[q], φ[q], mtable, lo, la, 0.25)
        p > 0 && mwet[p] && (depth[q] = D[p]; rh[q] = R[p]; global filled += 1)
    end
end
wet = depth .!= 0   # may be ≤ 0 (MOM6 wet down to -5 m); clipped at run time
depth[wet] .= max.(depth[wet], 1e-3)
@printf("MOM6 points %d (wet %d); our cells wet %d (ours %d, both %d, MOM6-only %d, ours-only %d), filled by nearest %d\n",
        length(D), count(mwet), count(wet), count(ourdepth .> 0), count(wet .& (ourdepth .> 0)), count(wet .& .!(ourdepth .> 0)),
        count(.!wet .& (ourdepth .> 0)), filled)
@printf("points per cell: mean %.2f, cells with 0: %d; rh max %.4f\n", sum(ntot) / length(ntot), count(ntot .== 0), maximum(rh))
jldsave(joinpath(@__DIR__, "..", "..", "GlobalKiloMischief", "tides", @sprintf("mom6_p04arc9_on_%d.jld2", Nx)); depth, rh)
