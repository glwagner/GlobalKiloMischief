# Offline comparison of self-attraction and loading potentials Ψ for one sea surface height snapshot of a 1/8ᵒ run,
# with pᵇ = g η over the ocean (the baroclinic part of the bottom pressure is not saved):
#
#   T40            degree 40 on a 48 × 96 Gaussian grid, evaluated exactly at cell centers (NumericalEarth's default)
#   T40 bilinear   the same coefficients, but Ψ on the Gaussian grid interpolated bilinearly to cell centers
#                  (closer to MITgcm, which interpolates bicubically)
#   T179           degree 179 on a 182 × 364 Gaussian grid (MITgcm's pkg/sal default resolution)
#   scalar         Ψ = -β g η with β = 0.1
#
#   julia --project=sal_env analysis/sal_offline.jl <output .jld2> <day>

using NumericalEarth
using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using SpeedyTransforms, ConservativeRegridding
using SpeedyTransforms.RingGrids: interpolate
using JLD2, Printf, Statistics, CairoMakie

file, day = ARGS[1], parse(Float64, ARGS[2])
g, β = 9.80665, 0.1

grid = TripolarGrid(CPU(); size = (2880, 1360, 1), halo = (7, 7, 4), z = (-6000, 0))
Nx, Ny, _ = size(grid)
bottom_height = regrid_bathymetry(grid; minimum_depth = 20, interpolation_passes = 3, major_basins = 1)
wet = interior(bottom_height, :, :, 1) .< 0

η, t = jldopen(file) do f
    iterations = keys(f["timeseries/t"])
    times = [f["timeseries/t/$i"] for i in iterations]
    n = argmin(abs.(times .- day * 86400))
    Float64.(f["timeseries/η/$(iterations[n])"][8:Nx+7, 8:Ny+7, 1]) .* wet, times[n] / 86400
end
@info @sprintf("η at day %.2f: max|η| = %.2f m", t, maximum(abs, η))

pᵇ = Field{Center, Center, Nothing}(grid)
interior(pᵇ, :, :, 1) .= g .* η

Ψ = Dict{String, Matrix{Float64}}()
for (name, kw) in (("T40", (; truncation = 40, nlat_half = 24)), ("T179", (; truncation = 179, nlat_half = 91)))
    build = @elapsed operator = NumericalEarth.Oceans.spherical_harmonic_loading_operator(grid; kw...)
    ψ = Field{Center, Center, Nothing}(grid)
    NumericalEarth.Oceans.apply_loading_operator!(ψ, operator, pᵇ)
    apply = @elapsed NumericalEarth.Oceans.apply_loading_operator!(ψ, operator, pᵇ)
    Ψ[name] = Array(interior(ψ, :, :, 1))
    @info @sprintf("%s: %s built in %.0f s, applied in %.2f s", name, summary(operator), build, apply)

    if name == "T40" # the weighted coefficients are still in operator.coefficients
        ψᵍ = SpeedyTransforms.transform(operator.coefficients, operator.transform)
        λ = mod.(vec(Array(λnodes(grid, Center(), Center()))), 360)
        φ = vec(Array(φnodes(grid, Center(), Center())))
        Ψ["T40 bilinear"] = - reshape(interpolate(λ, φ, ψᵍ), Nx, Ny)
    end
end
Ψ["scalar"] = - β * g .* η

# Gradients, which are what force the flow
function gradient_magnitude(ψ)
    f = Field{Center, Center, Nothing}(grid)
    interior(f, :, :, 1) .= ψ
    fill_halo_regions!(f)
    ∂xψ = Field(∂x(f))
    ∂yψ = Field(∂y(f))
    compute!(∂xψ); compute!(∂yψ)
    return sqrt.(interior(∂xψ, :, :, 1) .^ 2 .+ interior(∂yψ, :, :, 1) .^ 2)
end

rms(x) = sqrt(mean(abs2, x[wet]))
reference = Ψ["T40"]
∇reference = gradient_magnitude(reference)
for name in ("T40", "T40 bilinear", "T179", "scalar")
    δ = Ψ[name] .- reference
    ∇δ = gradient_magnitude(δ)
    @info @sprintf("%-13s rms Ψ = %.4f m² s⁻², rms(Ψ - Ψ_T40) / rms Ψ_T40 = %.3f, max|Ψ - Ψ_T40| = %.3f, rms|∇(Ψ - Ψ_T40)| / rms|∇Ψ_T40| = %.3f",
                   name, rms(Ψ[name]), rms(δ) / rms(reference), maximum(abs, δ[wet]), rms(∇δ) / rms(∇reference))
end

masked(x) = ifelse.(wet, x, NaN)
fig = Figure(size = (1600, 1100))
panels = (("Ψ, degree 40 (m² s⁻²)", reference),
          ("Ψ(scalar β = 0.1) − Ψ(degree 40)", Ψ["scalar"] .- reference),
          ("Ψ(degree 179) − Ψ(degree 40)", Ψ["T179"] .- reference),
          ("Ψ(bilinear) − Ψ(exact), degree 40", Ψ["T40 bilinear"] .- reference))
for (n, (label, field)) in enumerate(panels)
    ax = Axis(fig[(n - 1) ÷ 2 + 1, 2((n - 1) % 2) + 1]; title = label, aspect = DataAspect())
    limit = maximum(abs, field[wet])
    hm = heatmap!(ax, masked(field); colormap = :balance, colorrange = (-limit, limit))
    Colorbar(fig[(n - 1) ÷ 2 + 1, 2((n - 1) % 2) + 2], hm)
    hidedecorations!(ax)
end
Label(fig[0, :], @sprintf("Self-attraction and loading from η at day %.2f of the 1/8ᵒ run (index space; tripolar fold at top)", t))
filename = joinpath(@__DIR__, @sprintf("sal_offline_day%g.png", day))
save(filename, fig)
@info "Saved $filename"
