# Maps of the M2 tide from the harmonic analysis of the 1/8ᵒ runs (analysis/sal_harmonics.jld2, made by
# sal_harmonics.jl), compared with each other and with the EOT20 altimetric atlas (Hart-Davis et al. 2021).
#
#   julia --project=sal_env analysis/sal_m2_maps.jl
#
# Saves analysis/sal_m2_change.png (scalar SAL − no SAL and full SAL − no SAL) and analysis/sal_m2_vs_eot20.png
# (each run − EOT20), and prints area-weighted rms errors against EOT20.

using NumericalEarth
using Oceananigans
using Oceananigans.Grids: λnode, φnode
using Oceananigans.Operators: Azᶜᶜᶜ
using JLD2, CairoMakie, Printf, Statistics, Dates

const NCDatasets = Base.require(Base.PkgId(Base.UUID("85f8d34a-cbdd-5861-8df4-14fed0d494ab"), "NCDatasets"))

prefix = "internal_tides_CUDA_Float64_0.125deg_Nz10_dt3min_tau30d_"
runs = [(prefix * "catke", "no SAL"), (prefix * "scalar0.1_catke", "scalar SAL (β = 0.1)"),
        (prefix * "harmonic40_catke", "full SAL (degree 40)")]

results = jldopen(joinpath(@__DIR__, "sal_harmonics.jld2"))
amplitudes = results["amplitudes"]
first_day, last_day = results["first_day"], results["last_day"]
close(results)
Z(run, name) = amplitudes["$run/$name"]["M2"]

grid = TripolarGrid(CPU(); size = (2880, 1360, 1), halo = (7, 7, 4), z = (-6000, 0))
Nx, Ny = 2880, 1360
depth = -interior(regrid_bathymetry(grid; minimum_depth = 20, interpolation_passes = 3, major_basins = 1), :, :, 1)
λ = [mod(λnode(i, j, 1, grid, Center(), Center(), Center()), 360) for i in 1:Nx, j in 1:Ny]
φ = [φnode(i, j, 1, grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
area = [Azᶜᶜᶜ(i, j, 1, grid) for i in 1:Nx, j in 1:Ny]

#####
##### Model M2 as amplitude and Greenwich phase lag
#####

# The fit gives η ≈ Re(Z exp(iωt)) with t from 2015-01-01 00:00 UTC, and η = f A cos(ωt + V₀ + u - G), so
# A = |Z| / f and G = V₀ + u - arg Z. Astronomical arguments from mean longitudes (Meeus); u and f at the middle
# of the fitting window.
function julian_centuries(date)
    jd = Dates.datetime2julian(date)
    return (jd - 2451545.0) / 36525
end
T₀ = julian_centuries(DateTime(2015, 1, 1))
s = 218.3164477 + 481267.88123421T₀ # mean longitude of the Moon
h = 280.4664567 + 36000.76983T₀     # mean longitude of the Sun
V₀ = mod(2h - 2s, 360)              # M2 at 00:00 UT
Tₘ = julian_centuries(DateTime(2015, 1, 1) + Millisecond(round(Int, (first_day + last_day) / 2 * 86400e3)))
N = 125.04452 - 1934.136261Tₘ       # longitude of the Moon's ascending node
u = -2.14sind(N)
f = 1.0004 - 0.0373cosd(N) + 0.0002cosd(2N)
@info @sprintf("M2: V₀ = %.2fᵒ, u = %.2fᵒ, f = %.4f", V₀, u, f)

# Complex amplitude A exp(-iG) in metres
greenwich(Zᵐ) = @. abs(Zᵐ) / f * cis(-deg2rad(V₀ + u)) * cis(angle(Zᵐ))

#####
##### EOT20 at the model's cell centres (bilinear), amplitude and phase lag → A exp(-iG) in metres
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
observed = interpolate.(λ, φ) # NaN within one EOT20 cell of its coastline

#####
##### Errors against EOT20
#####

wet = depth .> 0
compared = wet .& isfinite.(observed)
deep = compared .& (depth .> 1000) .& (abs.(φ) .< 66)
area_rms(x, mask) = sqrt(sum(abs2.(x[mask]) .* area[mask]) / sum(area[mask]))
@info @sprintf("EOT20 rms M2 amplitude: %.2f cm (deep, |φ| < 66ᵒ), %.2f cm (all compared ocean)",
               100area_rms(abs.(observed), deep), 100area_rms(abs.(observed), compared))
modelled = Dict(label => greenwich(Z(run, "η")) for (run, label) in runs)
for (run, label) in runs
    m = modelled[label]
    # Best-fit uniform phase offset, to check the phase conventions: it should be small
    offset = rad2deg(angle(sum(m[deep] .* conj.(observed[deep]) .* area[deep])))
    @info @sprintf("%-22s rms M2 amplitude %.2f cm; vs EOT20, deep ocean: rms |ΔZ| %.2f cm, rms Δ amplitude %.2f cm, mean phase offset %.1fᵒ; all compared ocean: rms |ΔZ| %.2f cm",
                   label, 100area_rms(abs.(m), deep), 100area_rms(m .- observed, deep),
                   100area_rms(abs.(m) .- abs.(observed), deep), offset, 100area_rms(m .- observed, compared))
end

#####
##### Figures (grid index space; the grid is regular latitude-longitude south of about 55ᵒN)
#####

masked(x, mask) = ifelse.(mask, x, NaN)
function panel!(fig, row, col, field, mask, title; diverging = true, unit = "m", limit = quantile(abs.(field[mask]), 0.99))
    ax = Axis(fig[row, 2col - 1]; title, titlesize = 13, aspect = DataAspect())
    hm = heatmap!(ax, masked(field, mask); colormap = diverging ? :balance : :viridis,
                  colorrange = diverging ? (-limit, limit) : (0, limit))
    Colorbar(fig[row, 2col], hm; label = unit)
    hidedecorations!(ax)
end

reference_run = first(runs)[1]
wets = Dict(name => abs.(Z(reference_run, name)) .> 0 for name in ("η", "ζ"))
labels = Dict("η" => "sea surface height", "ζ" => "isopycnal displacement at 1177 m")

fig = Figure(size = (1700, 1250))
Label(fig[0, 1:4], @sprintf("M2 amplitude from a harmonic fit over days %g–%g", first_day, last_day); fontsize = 18)
for (col, name) in enumerate(("η", "ζ"))
    A₀ = abs.(Z(reference_run, name))
    panel!(fig, 1, col, A₀, wets[name], "no SAL: M2 $(labels[name])"; diverging = false)
    changes = [abs.(Z(run, name)) .- A₀ for (run, label) in runs[2:end]]
    limit = maximum(quantile(abs.(c[wets[name]]), 0.99) for c in changes) # one color scale for both rows
    for (row, (run, label)) in enumerate(runs[2:end])
        panel!(fig, row + 1, col, changes[row], wets[name], "$label − no SAL: M2 $(labels[name])"; limit)
    end
end
save(joinpath(@__DIR__, "sal_m2_change.png"), fig)
@info "Saved analysis/sal_m2_change.png"

fig = Figure(size = (1700, 1650))
Label(fig[0, 1:4], @sprintf("M2 sea surface height against EOT20 (fit over days %g–%g; ΔZ includes phase)", first_day, last_day); fontsize = 18)
panel!(fig, 1, 1, abs.(observed), compared, "EOT20 M2 amplitude"; diverging = false)
panel!(fig, 1, 2, rad2deg.(-angle.(observed)), compared, "EOT20 M2 Greenwich phase lag"; unit = "degrees")
# One color scale per column for the three runs, set by the run without SAL
amplitude_limit = quantile(abs.(abs.(modelled["no SAL"][compared]) .- abs.(observed[compared])), 0.99)
complex_limit = quantile(abs.(modelled["no SAL"][compared] .- observed[compared]), 0.99)
for (row, (run, label)) in enumerate(runs)
    m = modelled[label]
    panel!(fig, row + 1, 1, abs.(m) .- abs.(observed), compared, "$label − EOT20: amplitude"; limit = amplitude_limit)
    panel!(fig, row + 1, 2, abs.(m .- observed), compared, "|$label − EOT20|: complex difference"; diverging = false, limit = complex_limit)
end
save(joinpath(@__DIR__, "sal_m2_vs_eot20.png"), fig)
@info "Saved analysis/sal_m2_vs_eot20.png"
