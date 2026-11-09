# M2 from the single-layer runs (GlobalKiloMischief/tides/barotropic_tides.jl) with four treatments of
# self-attraction and loading, compared with the EOT20 altimetric atlas.
#
#   julia --project=sal_env analysis/barotropic_m2.jl [first day] [last day]
#
# Fits M2, S2, K1, O1 and a mean to hourly η over the window (default days 5–20), converts M2 to amplitude and
# Greenwich phase lag, and saves analysis/barotropic_m2.png: EOT20 on top, the four runs in the middle row, and
# each run − EOT20 below.

using NumericalEarth
using Oceananigans
using Oceananigans.Grids: λnode, φnode
using Oceananigans.Operators: Azᶜᶜᶜ
using JLD2, CairoMakie, Printf, Statistics, Dates, LinearAlgebra

const NCDatasets = Base.require(Base.PkgId(Base.UUID("85f8d34a-cbdd-5861-8df4-14fed0d494ab"), "NCDatasets"))

first_day = parse(Float64, get(ARGS, 1, "5"))
last_day = parse(Float64, get(ARGS, 2, "20"))
output = joinpath(@__DIR__, "..", "GlobalKiloMischief", "tides", "output")
# SUFFIX selects runs with options, e.g. SUFFIX=_wavedrag3
suffix = get(ENV, "SUFFIX", "")
# CELLS_PER_DEGREE (default 1 / 0.36) and DT (minutes, default 10) select the grid and time step
cells_per_degree = parse(Float64, get(ENV, "CELLS_PER_DEGREE", "2.7777777777777777"))
Δt_minutes = parse(Float64, get(ENV, "DT", "10"))
run_prefix = @sprintf("barotropic_tides_%gdeg_dt%gmin_", 1 / cells_per_degree, Δt_minutes)
filename(sal) = joinpath(output, "$run_prefix$sal$suffix.jld2")
runs = [("none", "no SAL"), ("scalar", "scalar SAL (β = 0.1)"), ("spherical_harmonic", "spherical harmonic SAL (n = 40)"),
        ("convolution_prem", "convolution SAL (PREM)")]
# Also compared with EOT20, but not drawn: the paper's convolution as in its MOM6 code (asymptotic Love numbers,
# no self-cell term, kernel ÷ 3), and as its equations are written (no ÷ 3)
extra_runs = [("convolution_paper", "paper's MOM6 convolution (÷3)"), ("convolution", "paper's equations (no ÷3)")]
extra_runs = filter(r -> isfile(filename(r[1])), extra_runs)
runs = filter(r -> isfile(filename(r[1])), runs)

Nx, Ny = round(Int, 360 * cells_per_degree), round(Int, 170 * cells_per_degree)
grid = TripolarGrid(CPU(); size = (Nx, Ny, 1), halo = (7, 7, 4), z = (-6000, 0))
depth = -interior(regrid_bathymetry(grid; minimum_depth = 20, interpolation_passes = 3, major_basins = 1), :, :, 1)
wet = depth .> 0
λ = [mod(λnode(i, j, 1, grid, Center(), Center(), Center()), 360) for i in 1:Nx, j in 1:Ny]
φ = [φnode(i, j, 1, grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
area = [Azᶜᶜᶜ(i, j, 1, grid) for i in 1:Nx, j in 1:Ny]

#####
##### Harmonic fit and Greenwich phase
#####

periods = (M2 = 12.4206012, S2 = 12.0, K1 = 23.9344696, O1 = 25.8193417) # hours
ω = [2π / (3600T) for T in values(periods)]

"""Complex M2 amplitude Z with η ≈ Re(Z exp(iωt)), t from 2015-01-01 00:00 UTC."""
function m2_fit(sal)
    jldopen(filename(sal)) do file
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
        @info @sprintf("%s: fitted %d hourly snapshots from day %.2f to %.2f", sal, length(window), τ[1] / 86400, τ[end] / 86400)
        return a .- im .* b
    end
end

julian_centuries(date) = (Dates.datetime2julian(date) - 2451545.0) / 36525
T₀ = julian_centuries(DateTime(2015, 1, 1))
V₀ = mod(2(280.4664567 + 36000.76983T₀) - 2(218.3164477 + 481267.88123421T₀), 360)
N = 125.04452 - 1934.136261julian_centuries(DateTime(2015, 1, 1) + Millisecond(round(Int, (first_day + last_day) / 2 * 86400e3)))
u = -2.14sind(N)
f = 1.0004 - 0.0373cosd(N) + 0.0002cosd(2N)

# A exp(-iG) in metres, with η = f A cos(ωt + V₀ + u - G)
greenwich(Z) = @. abs(Z) / f * cis(angle(Z) - deg2rad(V₀ + u))

modelled = Dict(sal => greenwich(m2_fit(sal)) .* wet for (sal, label) in [runs; extra_runs])

#####
##### EOT20 at cell centres
#####

eot = NCDatasets.NCDataset(joinpath(@__DIR__, "..", "egbert_ray", "data", "EOT20", "ocean_tides", "M2_ocean_eot20.nc"))
eot_λ = Float64.(eot["lon"][:]); eot_φ = Float64.(eot["lat"][:])
eot_Z = map((A, G) -> ismissing(A) || ismissing(G) ? complex(NaN) : A / 100 * cis(-deg2rad(G)), eot["amplitude"][:, :], eot["phase"][:, :])
close(eot)
Δλ, Δφ = eot_λ[2] - eot_λ[1], eot_φ[2] - eot_φ[1]
function interpolate(λ, φ)
    x = (λ - eot_λ[1]) / Δλ + 1; y = (φ - eot_φ[1]) / Δφ + 1
    i = clamp(floor(Int, x), 1, length(eot_λ) - 1); j = clamp(floor(Int, y), 1, length(eot_φ) - 1)
    a, b = x - i, y - j
    return (1 - a) * (1 - b) * eot_Z[i, j] + a * (1 - b) * eot_Z[i+1, j] + (1 - a) * b * eot_Z[i, j+1] + a * b * eot_Z[i+1, j+1]
end
observed = interpolate.(λ, φ)

compared = wet .& isfinite.(observed)
deep = compared .& (depth .> 1000) .& (abs.(φ) .< 66)
area_rms(x, mask) = sqrt(sum(abs2.(x[mask]) .* area[mask]) / sum(area[mask]))

errors = Dict{String, NamedTuple}()
@info @sprintf("EOT20 rms M2 amplitude %.2f cm (deep ocean: depth > 1000 m, |φ| < 66ᵒ)", 100area_rms(abs.(observed), deep))
for (sal, label) in [runs; extra_runs]
    m = modelled[sal]
    offset = rad2deg(angle(sum(m[deep] .* conj.(observed[deep]) .* area[deep])))
    errors[sal] = (deep = 100area_rms(m .- observed, deep), all = 100area_rms(m .- observed, compared),
                   amplitude = 100area_rms(abs.(m) .- abs.(observed), deep), offset = offset)
    @info @sprintf("%-32s rms amplitude %.2f cm; vs EOT20: rms |ΔZ| %.2f cm deep, %.2f cm all compared ocean; rms Δ amplitude %.2f cm deep; mean phase offset %.1fᵒ",
                   label, 100area_rms(abs.(m), deep), errors[sal].deep, errors[sal].all, errors[sal].amplitude, offset)
end

#####
##### Figure
#####

masked(x, mask) = ifelse.(mask, x, NaN)
fig = Figure(size = (2200, 1500))
Label(fig[0, 1:5], @sprintf("M2 amplitude in a single-layer model at %gᵒ%s (harmonic fit over days %g–%g), and EOT20", 1 / cells_per_degree, isempty(suffix) ? " without wave drag" : " with wave drag χ = " * replace(suffix, "_wavedrag" => ""), first_day, last_day);
      fontsize = 22)

amplitude_limit = quantile(abs.(observed[compared]), 0.99)
difference_limit = quantile(abs.(abs.(modelled[first(runs)[1]][compared]) .- abs.(observed[compared])), 0.99)

top = fig[1, 2:3] = GridLayout()
ax = Axis(top[1, 1]; title = "EOT20 (observed)", titlesize = 18, aspect = DataAspect())
hm = heatmap!(ax, masked(abs.(observed), compared); colormap = :viridis, colorrange = (0, amplitude_limit))
Colorbar(top[1, 2], hm; label = "M2 amplitude (m)")
hidedecorations!(ax)

for (col, (sal, label)) in enumerate(runs)
    ax = Axis(fig[2, col]; title = label, titlesize = 15, aspect = DataAspect())
    hm = heatmap!(ax, masked(abs.(modelled[sal]), compared); colormap = :viridis, colorrange = (0, amplitude_limit))
    col == 4 && Colorbar(fig[2, 5], hm; label = "M2 amplitude (m)")
    hidedecorations!(ax)

    e = errors[sal]
    ax = Axis(fig[3, col]; title = @sprintf("%s − EOT20\nrms error %.1f cm (deep ocean, with phase)", label, e.deep),
              titlesize = 15, aspect = DataAspect())
    hm = heatmap!(ax, masked(abs.(modelled[sal]) .- abs.(observed), compared); colormap = :balance,
                  colorrange = (-difference_limit, difference_limit))
    col == 4 && Colorbar(fig[3, 5], hm; label = "amplitude difference (m)")
    hidedecorations!(ax)
end
rowsize!(fig.layout, 1, Relative(0.42))
figure_name = @sprintf("barotropic_m2_%gdeg_dt%gmin%s.png", 1 / cells_per_degree, Δt_minutes, suffix)
save(joinpath(@__DIR__, figure_name), fig)
@info "Saved analysis/$figure_name"
