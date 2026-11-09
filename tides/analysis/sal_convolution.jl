# Reproduce Figs. 1–3 of Chen, Wang, Arbic & Krasny (2026, arXiv:2602.01416) at 0.36ᵒ, and separate the causes of the
# differences between spherical harmonic and convolution self-attraction and loading (SAL).
#
#   julia -t auto --project=sal_env analysis/sal_convolution.jl <1/8ᵒ run> <day>
#
# η from one snapshot of a 1/8ᵒ run is averaged conservatively onto a 1000 × 472 TripolarGrid (0.36ᵒ), with
# minimum depth 20 m. The SAL acceleration g∇η_SAL = -∇Ψ is computed with
#
#   SH n          spherical harmonics to degree n by quadrature over the ocean cells, as MOM6 does
#                 ("exact": Gegout load Love numbers; "asymptotic": the paper's fit k′ = a₁/n, h′ = b₀ + b₁/n)
#   convolution   NumericalEarth's ConvolutionLoadingOperator (asymptotic Love numbers, i = j term dropped)
#   hybrid        convolution + SH to degree 40 of (exact - asymptotic) Love numbers
#
# at cell centres, eastward and northward, with the SH gradient evaluated analytically. The paper's discretization
# (SH potential differenced to faces, convolution acceleration averaged to faces) is computed as well. Saves
# analysis/sal_convolution_fig1.png, _fig2.png, _fig3.png, _diagnostics.png and _results.jld2.

using NumericalEarth
using Oceananigans
using Oceananigans.Grids: λnode, φnode
using Oceananigans.Operators: Azᶜᶜᶜ, Δxᶠᶜᶜ, Δyᶜᶠᶜ
using Oceananigans.ImmersedBoundaries: GridFittedBottom
using SpeedyTransforms, ConservativeRegridding
using ConservativeRegridding: Regridder, regrid!
using JLD2, CairoMakie, Printf, Statistics, LinearAlgebra
using Base.Threads: @threads, nthreads, threadid

const Ext = Base.get_extension(NumericalEarth, :NumericalEarthSpeedyTransformsExt)

run_name, day = ARGS[1], parse(Float64, ARGS[2])
output = joinpath(@__DIR__, "..", "GlobalKiloMischief", "tides", "output")
g = 9.80665
ρ₀, ρₑ = 1035, 5517
@info "$(nthreads()) threads"

#####
##### η on the 0.36ᵒ grid
#####

source_grid = TripolarGrid(CPU(); size = (2880, 1360, 1), halo = (7, 7, 4), z = (-6000, 0))
grid_with_depth = TripolarGrid(CPU(); size = (1000, 472, 1), halo = (7, 7, 4), z = (-6000, 0))
source_wet = interior(regrid_bathymetry(source_grid; minimum_depth = 20, interpolation_passes = 3, major_basins = 1), :, :, 1) .< 0
bottom_height = regrid_bathymetry(grid_with_depth; minimum_depth = 20, interpolation_passes = 3, major_basins = 1)
Nx, Ny = 1000, 472
wet = interior(bottom_height, :, :, 1) .< 0

η₈, t = jldopen(joinpath(output, run_name * ".jld2")) do file
    iterations = keys(file["timeseries/t"])
    times = [file["timeseries/t/$i"] for i in iterations]
    n = argmin(abs.(times .- day * 86400))
    Float64.(file["timeseries/η/$(iterations[n])"][8:2887, 8:1367, 1]) .* source_wet, times[n] / 86400
end

regrid_time = @elapsed begin
    destination = Ext.horizontal_grid(grid_with_depth)
    regridder = Regridder(ConservativeRegridding.GOCore.best_manifold(destination), destination, Ext.horizontal_grid(source_grid))
    fraction = zeros(Nx * Ny); weighted = zeros(Nx * Ny)
    regrid!(fraction, regridder, vec(Float64.(source_wet)))
    regrid!(weighted, regridder, vec(η₈))
end
fraction = reshape(fraction, Nx, Ny); weighted = reshape(weighted, Nx, Ny)
η = @. ifelse(wet & (fraction > 1e-3), weighted / fraction, 0.0)
@info @sprintf("η at day %.2f regridded to 0.36ᵒ in %.0f s: %d wet cells, max|η| = %.2f m, %d wet cells without source ocean",
               t, regrid_time, count(wet), maximum(abs, η), count(wet .& (fraction .≤ 1e-3)))

# A grid whose only active cells are the wet columns, for the convolution operator
grid = ImmersedBoundaryGrid(TripolarGrid(CPU(); size = (Nx, Ny, 1), halo = (7, 7, 4), z = (-1, 0)),
                            GridFittedBottom(ifelse.(wet, -10.0, 0.0)))
R = grid.radius
λᶜ = [λnode(i, j, 1, grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
φᶜ = [φnode(i, j, 1, grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
Ω = [Azᶜᶜᶜ(i, j, 1, grid) / R^2 for i in 1:Nx, j in 1:Ny]
pᵇ = g .* η

#####
##### Load Love numbers
#####

exact_love = Ext.load_love_numbers(nothing)
asymptotic = AsymptoticLoadLoveNumbers()
γ_exact(n) = 3ρ₀ / ρₑ * (1 + exact_love.k[n+1] - exact_love.h[n+1]) / (2n + 1)
γ_asymptotic(n) = 3ρ₀ / ρₑ * NumericalEarth.Oceans.one_plus_k_minus_h(asymptotic, n) / (2n + 1)

for n in (0, 1, 2, 3, 4, 5, 10, 20, 40, 100, 200, 500, 1000)
    exact = 1 + exact_love.k[n+1] - exact_love.h[n+1]
    fit = NumericalEarth.Oceans.one_plus_k_minus_h(asymptotic, n)
    @info @sprintf("n = %4d: 1 + k′ - h′ = %.4f (Gegout, CM), %.4f (asymptotic fit), ratio %.3f", n, exact, fit, fit / exact)
end

#####
##### Spherical harmonics by quadrature on arbitrary points (orthonormal on the unit sphere, stored by m then ℓ)
#####

index_offsets(L) = cumsum([0; [L + 1 - m for m in 0:L-1]])

"""Coefficients cₗₘ = Σ f Ω conj(Yₗₘ) over the points where mask is true."""
function analyze(f, λ, φ, Ω, mask, L)
    recursion = Ext.legendre_recursion_coefficients(Float64, L)
    points = findall(mask)
    partial = [zeros(ComplexF64, (L + 1) * (L + 2) ÷ 2) for _ in 1:Threads.maxthreadid()]
    chunks = collect(Iterators.partition(points, cld(length(points), 4nthreads())))
    @threads :static for chunk in chunks
        c = partial[threadid()]
        for p in chunk
            x, y = sind(φ[p]), cosd(φ[p])
            w = f[p] * Ω[p]
            e = cis(-deg2rad(λ[p]))
            eᵐ = complex(w)
            Pᵐᵐ = 1.0
            n = 1
            for m in 0:L
                Pᵐᵐ *= recursion[n, 1] * (m == 0 ? 1.0 : y)
                P², P¹ = 0.0, Pᵐᵐ
                c[n] += eᵐ * Pᵐᵐ
                for q in n+1:n+L-m
                    P = recursion[q, 1] * (x * P¹ - recursion[q, 2] * P²)
                    c[q] += eᵐ * P
                    P², P¹ = P¹, P
                end
                eᵐ *= e
                n += L + 1 - m
            end
        end
    end
    return sum(partial)
end

# Σₘ Sₘ(φ) and its λ derivative, Sₘ = Σₗ cₗₘ P̄ₗₘ(sin φ)
function synthesize_point(c, recursion, L, λ, φ)
    x, y = sind(φ), cosd(φ)
    ψ = ∂λψ = 0.0
    Pᵐᵐ = 1.0
    n = 1
    e = cis(deg2rad(λ)); eᵐ = complex(1.0)
    for m in 0:L
        Pᵐᵐ *= recursion[n, 1] * (m == 0 ? 1.0 : y)
        P², P¹ = 0.0, Pᵐᵐ
        S = c[n] * Pᵐᵐ
        for q in n+1:n+L-m
            P = recursion[q, 1] * (x * P¹ - recursion[q, 2] * P²)
            S += c[q] * P
            P², P¹ = P¹, P
        end
        factor = m == 0 ? 1 : 2
        ψ += factor * real(S * eᵐ)
        ∂λψ += factor * real(im * m * S * eᵐ)
        eᵐ *= e
        n += L + 1 - m
    end
    return ψ, ∂λψ
end

"""Potential and its eastward and northward gradient (per metre) of Σ cₗₘ Yₗₘ at cell centres."""
function synthesize(c, L, λ, φ; h = 1e-3)
    recursion = Ext.legendre_recursion_coefficients(Float64, L)
    ψ = zeros(size(λ)); ∂x = zeros(size(λ)); ∂y = zeros(size(λ))
    @threads for p in eachindex(λ)
        ψ[p], ∂λψ = synthesize_point(c, recursion, L, λ[p], φ[p])
        ψ⁺, _ = synthesize_point(c, recursion, L, λ[p], φ[p] + h)
        ψ⁻, _ = synthesize_point(c, recursion, L, λ[p], φ[p] - h)
        ∂x[p] = ∂λψ / (R * cosd(φ[p]))
        ∂y[p] = (ψ⁺ - ψ⁻) / (2deg2rad(h) * R)
    end
    return ψ, ∂x, ∂y
end

function weight!(c, L, γ)
    n = 1
    for m in 0:L, ℓ in m:L
        c[n] *= γ(ℓ)
        n += 1
    end
    return c
end

truncate_coefficients(c, L, L′) = [c[index_offsets(L)[m+1] + ℓ - m + 1] for m in 0:L′ for ℓ in m:L′]

"""Power per degree of the coefficients."""
function degree_power(c, L)
    power = zeros(L + 1)
    n = 1
    for m in 0:L, ℓ in m:L
        power[ℓ+1] += (m == 0 ? 1 : 2) * abs2(c[n])
        n += 1
    end
    return power
end

#####
##### SAL accelerations at cell centres: a = -∇Ψ = ∇(Σ γₗ pᵇₗₘ Yₗₘ)
#####

Lmax = 400
analysis_time = @elapsed p̂ = analyze(pᵇ, λᶜ, φᶜ, Ω, wet, Lmax)
@info @sprintf("pᵇ analysed to degree %d in %.0f s", Lmax, analysis_time)

accelerations = Dict{String, NTuple{3, Matrix{Float64}}}() # (Ψ, eastward, northward)
for (name, L, γ) in (("SH 40", 40, γ_exact), ("SH 200", 200, γ_exact),
                     ("SH 40 asymptotic", 40, γ_asymptotic), ("SH 200 asymptotic", 200, γ_asymptotic),
                     ("SH 400 asymptotic", 400, γ_asymptotic),
                     ("hybrid correction", 40, n -> γ_exact(n) - γ_asymptotic(n)))
    c = weight!(truncate_coefficients(p̂, Lmax, L), L, γ)
    time = @elapsed ψ, ax, ay = synthesize(c, L, λᶜ, φᶜ)
    accelerations[name] = (-ψ, ax, ay)
    @info @sprintf("%s synthesized in %.0f s", name, time)
end

operator = ConvolutionLoadingOperator(grid)
p_field = Field{Center, Center, Nothing}(grid); interior(p_field, :, :, 1) .= pᵇ
∂xΨ = Field{Center, Center, Nothing}(grid); ∂yΨ = Field{Center, Center, Nothing}(grid)
convolution_time = @elapsed loading_gradient!(∂xΨ, ∂yΨ, operator, p_field)
@info @sprintf("%s applied in %.0f s", summary(operator), convolution_time)
convolution = (zeros(Nx, Ny), -Array(interior(∂xΨ, :, :, 1)) .* wet, -Array(interior(∂yΨ, :, :, 1)) .* wet)
accelerations["convolution"] = convolution
correction = accelerations["hybrid correction"]
accelerations["hybrid"] = (correction[1], convolution[2] .+ correction[2], convolution[3] .+ correction[3])

#####
##### The paper's discretization: Ψ differenced to faces (SH), centre accelerations averaged to faces (convolution)
#####

# Local angle of the grid's i direction from east, from neighbouring centres
θ = [begin
         i⁺, i⁻ = mod1(i + 1, Nx), mod1(i - 1, Nx)
         dλ = mod(λᶜ[i⁺, j] - λᶜ[i⁻, j] + 180, 360) - 180
         atan(φᶜ[i⁺, j] - φᶜ[i⁻, j], dλ * cosd(φᶜ[i, j]))
     end for i in 1:Nx, j in 1:Ny]

u_face_wet = [wet[i, j] && wet[mod1(i - 1, Nx), j] for i in 1:Nx, j in 1:Ny]
v_face_wet = [j > 1 && wet[i, j] && wet[i, j-1] for i in 1:Nx, j in 1:Ny]
Δx = [Δxᶠᶜᶜ(i, j, 1, grid) for i in 1:Nx, j in 1:Ny]
Δy = [Δyᶜᶠᶜ(i, j, 1, grid) for i in 1:Nx, j in 1:Ny]

function differenced_to_faces(Ψ)
    au = [u_face_wet[i, j] ? -(Ψ[i, j] - Ψ[mod1(i - 1, Nx), j]) / Δx[i, j] : 0.0 for i in 1:Nx, j in 1:Ny]
    av = [v_face_wet[i, j] ? -(Ψ[i, j] - Ψ[i, j-1]) / Δy[i, j] : 0.0 for i in 1:Nx, j in 1:Ny]
    return au, av
end

function averaged_to_faces(ae, an)
    ax = @. ae * cos(θ) + an * sin(θ) # grid-aligned components
    ay = @. -ae * sin(θ) + an * cos(θ)
    au = [u_face_wet[i, j] ? (ax[i, j] + ax[mod1(i - 1, Nx), j]) / 2 : 0.0 for i in 1:Nx, j in 1:Ny]
    av = [v_face_wet[i, j] ? (ay[i, j] + ay[i, j-1]) / 2 : 0.0 for i in 1:Nx, j in 1:Ny]
    return au, av
end

faces = Dict{String, NTuple{2, Matrix{Float64}}}()
for name in ("SH 40", "SH 200")
    faces[name] = differenced_to_faces(accelerations[name][1])
end
faces["convolution"] = averaged_to_faces(convolution[2], convolution[3])
faces["SH 200, averaged"] = averaged_to_faces(accelerations["SH 200"][2], accelerations["SH 200"][3])

#####
##### Diagnostics
#####

# Distance from land in cells (8-connected), for separating coastal from open-ocean differences
distance = fill(typemax(Int), Nx, Ny)
distance[.!wet] .= 0
for d in 1:12
    for j in 1:Ny, i in 1:Nx
        distance[i, j] == typemax(Int) || continue
        near = any(distance[mod1(i + di, Nx), clamp(j + dj, 1, Ny)] == d - 1 for di in -1:1, dj in -1:1)
        near && (distance[i, j] = d)
    end
end

open_ocean = wet .& (distance .> 3) .& (abs.(φᶜ) .< 66)
area_rms(a, mask) = sqrt(sum(abs2.(a[mask]) .* Ω[mask]) / sum(Ω[mask]))
magnitude(name) = sqrt.(accelerations[name][2] .^ 2 .+ accelerations[name][3] .^ 2)

reference = "SH 400 asymptotic"
for name in ("SH 40", "SH 200", "SH 40 asymptotic", "SH 200 asymptotic", "SH 400 asymptotic", "convolution", "hybrid")
    a = accelerations[name]; r = accelerations[reference]
    difference = sqrt.((a[2] .- r[2]) .^ 2 .+ (a[3] .- r[3]) .^ 2)
    @info @sprintf("%-18s rms |a|: all ocean %.3e, within 3 cells of land %.3e, open ocean %.3e m s⁻²; rms |a - a(%s)|: all %.3e, open ocean %.3e",
                   name, area_rms(magnitude(name), wet), area_rms(magnitude(name), wet .& (distance .≤ 3)),
                   area_rms(magnitude(name), open_ocean), reference, area_rms(difference, wet), area_rms(difference, open_ocean))
end
for name in keys(faces)
    au, av = faces[name]
    @info @sprintf("faces, %-17s rms u-face acceleration %.3e, v-face %.3e m s⁻²", name, area_rms(au, u_face_wet), area_rms(av, v_face_wet))
end

#####
##### Spectra of the zonal acceleration over ocean points (Fig. 3), unnormalized and normalized
#####

Lspectrum = 400
spectra = Dict{String, Vector{Float64}}()
for name in ("SH 40", "SH 200", "SH 400 asymptotic", "convolution", "hybrid")
    spectra[name] = degree_power(analyze(accelerations[name][2], λᶜ, φᶜ, Ω, wet, Lspectrum), Lspectrum)
    spectra[name * ", open ocean only"] = degree_power(analyze(accelerations[name][2], λᶜ, φᶜ, Ω, open_ocean, Lspectrum), Lspectrum)
end
λᵘ = [λnode(i, j, 1, grid, Face(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
φᵘ = [φnode(i, j, 1, grid, Face(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
for name in ("SH 40", "SH 200", "convolution")
    spectra[name * ", faces"] = degree_power(analyze(faces[name][1], λᵘ, φᵘ, Ω, u_face_wet, Lspectrum), Lspectrum)
end
@info "Spectra done"

for (band, degrees) in (("n = 2–10", 3:11), ("n = 11–40", 12:41), ("n = 41–200", 42:201), ("n = 201–400", 202:401))
    line = join([@sprintf("%s %.2e", name, sum(spectra[name][degrees])) for name in
                 ("SH 40", "SH 200", "SH 400 asymptotic", "convolution", "hybrid", "SH 200, faces", "convolution, faces",
                  "SH 200, open ocean only", "convolution, open ocean only")], "; ")
    @info "zonal acceleration power, $band: $line"
end

jldsave(joinpath(@__DIR__, "sal_convolution_results.jld2"); run_name, t, η, wet, λᶜ, φᶜ, Ω, distance,
        accelerations, faces, spectra, p̂, Lmax)

#####
##### Figures
#####

masked(a, mask) = ifelse.(mask, a, NaN)

function map_figure(names, field_sets, filename, title; is = 1:Nx, js = 1:Ny)
    fig = Figure(size = (1600, 420 * length(names) + 80))
    Label(fig[0, 1:4], title; fontsize = 18)
    limit = 0.0
    for (r, (name, (au, av, mu, mv))) in enumerate(zip(names, field_sets)), (col, (a, m, component)) in enumerate(((au, mu, "zonal"), (av, mv, "meridional")))
        r == 1 && col == 1 && (limit = quantile(abs.(a[is, js][m[is, js]]), 0.995))
        ax = Axis(fig[r, 2col - 1]; title = "$name, $component", aspect = DataAspect())
        hm = heatmap!(ax, is, js, masked(a, m)[is, js]; colormap = :balance, colorrange = (-limit, limit))
        Colorbar(fig[r, 2col], hm; label = "m s⁻²")
        hidedecorations!(ax)
    end
    save(filename, fig)
    @info "Saved $filename"
end

panel_names = ("SH n = 40 (differenced to faces)", "SH n = 200 (differenced to faces)", "Convolution (averaged to faces)")
sets = [(faces[n]..., u_face_wet, v_face_wet) for n in ("SH 40", "SH 200", "convolution")]
day_title = @sprintf("g∇η_SAL from η at day %.2f of %s, regridded to 0.36ᵒ (grid index space; regular latitude-longitude south of 55ᵒN)", t, run_name)
map_figure(panel_names, sets, joinpath(@__DIR__, "sal_convolution_fig1.png"), day_title)

# Fig. 2: Northeast Atlantic and North Sea
box = findall(@. (-20 ≤ mod(λᶜ + 180, 360) - 180 ≤ 15) & (45 ≤ φᶜ ≤ 65))
is = minimum(p[1] for p in box):maximum(p[1] for p in box)
js = minimum(p[2] for p in box):maximum(p[2] for p in box)
map_figure(panel_names, sets, joinpath(@__DIR__, "sal_convolution_fig2.png"), "As Fig. 1, Northeast Atlantic and North Sea"; is, js)

# Centre comparison with the same Love numbers, to separate Gibbs ringing from the discretization
centre_names = ("SH n = 200, asymptotic Love numbers, centres", "Convolution, centres", "Convolution − SH n = 400 (asymptotic), centres")
r = accelerations["SH 400 asymptotic"]
centre_sets = [(accelerations["SH 200 asymptotic"][2], accelerations["SH 200 asymptotic"][3], wet, wet),
               (convolution[2], convolution[3], wet, wet),
               (convolution[2] .- r[2], convolution[3] .- r[3], wet, wet)]
map_figure(centre_names, centre_sets, joinpath(@__DIR__, "sal_convolution_centres.png"),
           "Same Love numbers, all at cell centres (eastward and northward)")

fig = Figure(size = (1500, 600))
n = 1:Lspectrum
for (col, (normalize, label)) in enumerate(((false, "power (m² s⁻⁴)"), (true, "power, normalized to unit sum")))
    ax = Axis(fig[1, col]; xscale = log10, yscale = log10, xlabel = "spherical harmonic degree n", ylabel = label,
              title = "Zonal SAL acceleration over ocean points" * (normalize ? " (as the paper's Fig. 3)" : ""))
    for (name, style) in (("SH 40", (color = :steelblue,)), ("SH 200", (color = :darkorange,)),
                          ("convolution", (color = :black,)), ("hybrid", (color = :purple, linestyle = :dot)),
                          ("SH 400 asymptotic", (color = :green, linestyle = :dash)),
                          ("SH 200, faces", (color = :darkorange, linestyle = :dashdot)),
                          ("convolution, faces", (color = :gray, linestyle = :dashdot)),
                          ("SH 200, open ocean only", (color = :darkorange, linewidth = 0.8, linestyle = :dot)),
                          ("convolution, open ocean only", (color = :black, linewidth = 0.8, linestyle = :dot)))
        s = spectra[name][2:end]
        normalize && (s = s ./ sum(spectra[name]))
        lines!(ax, n, max.(s, 1e-30); label = name, style...)
    end
    col == 2 && axislegend(ax; position = :lb, labelsize = 10)
end
save(joinpath(@__DIR__, "sal_convolution_fig3.png"), fig)
@info "Saved analysis/sal_convolution_fig3.png"
