# M2 vs EOT20 for any set of single-layer runs, fitting M2, N2, S2, K1, O1 and a mean (N2 separates from M2 only in
# windows ≥ 27.6 d; barotropic_m2.jl's 15-day fit without N2 lets N2 leak into M2).
#   julia --project=sal_env analysis/diag/m2_table.jl <first day> <last day> <figure name or -> <run name>...
# Run names are output files without .jld2. CELLS_PER_DEGREE as in barotropic_m2.jl. FIT_N2=false reproduces the old fit.
using NumericalEarth, Oceananigans
using Oceananigans.Grids: λnode, φnode
using Oceananigans.Operators: Azᶜᶜᶜ
using JLD2, CairoMakie, Printf, Statistics, Dates, LinearAlgebra
const NCDatasets = Base.require(Base.PkgId(Base.UUID("85f8d34a-cbdd-5861-8df4-14fed0d494ab"), "NCDatasets"))

first_day, last_day = parse(Float64, ARGS[1]), parse(Float64, ARGS[2])
figure_name = ARGS[3]
names = ARGS[4:end]
output = joinpath(@__DIR__, "..", "..", "GlobalKiloMischief", "tides", "output")
cells_per_degree = parse(Float64, get(ENV, "CELLS_PER_DEGREE", "2.7777777777777777"))
Nx, Ny = round(Int, 360 * cells_per_degree), round(Int, 170 * cells_per_degree)
grid = TripolarGrid(CPU(); size = (Nx, Ny, 1), halo = (7, 7, 4), z = (-6000, 0))
depth = -interior(regrid_bathymetry(grid; minimum_depth = 20, interpolation_passes = 3, major_basins = 1), :, :, 1)
wet = depth .> 0
λ = [mod(λnode(i, j, 1, grid, Center(), Center(), Center()), 360) for i in 1:Nx, j in 1:Ny]
φ = [φnode(i, j, 1, grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
area = [Azᶜᶜᶜ(i, j, 1, grid) for i in 1:Nx, j in 1:Ny]

periods = get(ENV, "FIT_N2", "true") == "true" ? [12.4206012, 12.65834751, 12.0, 23.9344696, 25.8193417] :
                                                 [12.4206012, 12.0, 23.9344696, 25.8193417]
ω = [2π / (3600T) for T in periods]
function m2_fit(name)
    if endswith(name, "_harmonics")
        h = load(joinpath(output, name * ".jld2"))
        @info @sprintf("%s: online fit of %d samples, days %.2f–%.2f", name, h["samples"], h["first_day"], h["last_day"])
        global window_mid = (h["first_day"] + h["last_day"]) / 2
        return h["Z"]["M2"]
    end
    global window_mid = (first_day + last_day) / 2
    jldopen(joinpath(output, name * ".jld2")) do file
        iterations = keys(file["timeseries/t"])
        times = [file["timeseries/t/$i"] for i in iterations]
        window = findall(t -> first_day * 86400 - 1 ≤ t ≤ last_day * 86400 + 1, times)
        τ = times[window]
        A = hcat(ones(length(τ)), [cos(w * t) for t in τ, w in ω], [sin(w * t) for t in τ, w in ω])
        P = pinv(A); K = length(ω)
        a = zeros(Nx, Ny); b = zeros(Nx, Ny)
        for (n, i) in enumerate(iterations[window])
            η = Float64.(file["timeseries/η/$i"][8:Nx+7, 8:Ny+7, 1])
            a .+= P[2, n] .* η; b .+= P[2 + K, n] .* η
        end
        @info @sprintf("%s: %d snapshots, days %.2f–%.2f", name, length(window), τ[1] / 86400, τ[end] / 86400)
        return a .- im .* b
    end
end
julian_centuries(date) = (Dates.datetime2julian(date) - 2451545.0) / 36525
T₀ = julian_centuries(DateTime(2015, 1, 1))
V₀ = mod(2(280.4664567 + 36000.76983T₀) - 2(218.3164477 + 481267.88123421T₀), 360)
window_mid = (first_day + last_day) / 2
function greenwich(Z)
    N = 125.04452 - 1934.136261julian_centuries(DateTime(2015, 1, 1) + Millisecond(round(Int, window_mid * 86400e3)))
    u = -2.14sind(N); f = 1.0004 - 0.0373cosd(N) + 0.0002cosd(2N)
    return @. abs(Z) / f * cis(angle(Z) - deg2rad(V₀ + u))
end

eot = NCDatasets.NCDataset(joinpath(@__DIR__, "..", "..", "egbert_ray", "data", "EOT20", "ocean_tides", "M2_ocean_eot20.nc"))
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
rms(x, mask) = sqrt(sum(abs2.(x[mask]) .* area[mask]) / sum(area[mask]))
# Total rms of |Z_m - Z_o| (= rms of η difference × √2); amplitude part rms(|Z_m| - |Z_o|); phase part the rest
sphere = 4π * 6371e3^2
fdeep = sum(area[deep]) / sphere; fall = sum(area[compared]) / sphere
println(@sprintf("deep area = %.3f, compared area = %.3f of the sphere", fdeep, fall))
println("| run (days $(first_day)–$(last_day)) | deep total | deep amp | deep phase | global total | deep total, ÷4πR² (paper eq. 22 literally) | global ÷4πR² | rms amp deep | phase offset (ᵒ, + = model leads) |")
println("|---|---|---|---|---|---|---|---|---|")
println(@sprintf("| EOT20 | | | | | | | %.2f | |", 100rms(abs.(observed), deep)))
models = Dict{String, Matrix{ComplexF64}}()
for name in names
    Zfit = m2_fit(name); m = greenwich(Zfit) .* wet
    models[name] = m
    # The paper's D² = ½|ΔZ|², so its errors are rms |ΔZ| / √2 (cm)
    tot = 100rms(m .- observed, deep) / √2; amp = 100rms(abs.(m) .- abs.(observed), deep) / √2
    glob = 100rms(m .- observed, compared) / √2
    off = rad2deg(angle(sum(m[deep] .* conj.(observed[deep]) .* area[deep])))
    println(@sprintf("| %s | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f | %.1f |", replace(name, "barotropic_tides_0.36deg_" => "", "_harmonics" => ""),
                     tot, amp, sqrt(max(tot^2 - amp^2, 0)), glob, tot * sqrt(fdeep), glob * sqrt(fall), 100rms(abs.(m), deep), off))
end
flush(stdout)
if figure_name != "-"
    masked(x, mask) = ifelse.(mask, x, NaN)
    n = length(names)
    fig = Figure(size = (700n, 900))
    lim = quantile(abs.(observed[compared]), 0.99)
    for (c, name) in enumerate(names)
        m = models[name]
        ax = Axis(fig[1, c]; title = name, titlesize = 11, aspect = DataAspect()); hidedecorations!(ax)
        hm = heatmap!(ax, masked(abs.(m), compared); colormap = :viridis, colorrange = (0, lim))
        c == n && Colorbar(fig[1, n+1], hm; label = "M2 amplitude (m)")
        ax = Axis(fig[2, c]; title = @sprintf("|Z − Z_EOT20|/√2, deep rms %.1f cm", 100rms(m .- observed, deep) / √2), titlesize = 11, aspect = DataAspect()); hidedecorations!(ax)
        hm = heatmap!(ax, masked(abs.(m .- observed) ./ √2, compared); colormap = :amp, colorrange = (0, 0.3))
        c == n && Colorbar(fig[2, n+1], hm; label = "error (m)")
        ax = Axis(fig[3, c]; title = "phase − EOT20 phase (ᵒ, + = model leads)", titlesize = 11, aspect = DataAspect()); hidedecorations!(ax)
        hm = heatmap!(ax, masked(rad2deg.(angle.(m .* conj.(observed))), compared .& (abs.(observed) .> 0.05)); colormap = :balance, colorrange = (-90, 90))
        c == n && Colorbar(fig[3, n+1], hm; label = "ᵒ")
    end
    save(joinpath(@__DIR__, figure_name), fig)
    @info "saved analysis/diag/$figure_name"
end
