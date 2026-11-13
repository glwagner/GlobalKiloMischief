# Two globe views of surface ocean speed with sea ice speed on top, from the per-rank surface and sea ice output of a
# distributed run: shared by the figure (visualize_global_simulation_twelfth.jl) and the animation
# (animate_global_simulation_twelfth.jl). Adapted from OceananigansPaper figures/visualize_global_simulation3.jl:
#
# - Speed from face values averaged to cell centers: s = √((uᵢ² + uᵢ₊₁²)/2 + (vⱼ² + vⱼ₊₁²)/2). The original summed the
#   two faces, which made speeds √2 too large.
# - Sea ice is drawn where its concentration ℵ ≥ 1e-3; elsewhere the ocean shows through.
# - The surface mesh is closed across the periodic longitude seam and across the north fold.
# - The grid starts at 80ᵒS: south of it is a plain neutral gray cap, with no data and no label.
# - Fixed physical color ranges instead of fractions of the maxima.

using CairoMakie
using JLD2
using Oceananigans
using Dates
using Printf

const start_date = DateTime(2015, 1, 1)
const ocean_speed_range = (0, 1.0)   # m s⁻¹
const ice_speed_range = (0, 0.4)     # m s⁻¹
const minimum_ice_concentration = 1e-3

rank_files(directory, prefix) = [joinpath(directory, f) for f in sort(filter(f -> startswith(f, prefix * "_rank"), readdir(directory)),
                                                                      by = f -> parse(Int, match(r"rank(\d+)", f)[1]))]

# Iterations and times of a ranked output, in time order
function snapshots(directory, prefix)
    jldopen(first(rank_files(directory, prefix))) do file
        iterations = sort(filter(!=("serialized"), collect(keys(file["timeseries/t"]))), by = x -> parse(Int, x))
        [(iteration, file["timeseries/t/$iteration"]) for iteration in iterations]
    end
end

# Stack the interiors of the ranks' bands of rows. Each array has a halo of `(size % 360) / 2` cells in x (the grid is
# 360c × 170c); the y halo can be wider than 170 rows (the sea ice carries 2 × substeps + 3), so the band height is
# taken from the global Ny = 170c instead. Faces keep one more column or row than centers.
function stitch(directory, prefix, name, iteration; face_x = false, face_y = false)
    paths = rank_files(directory, prefix)
    bands = map(paths) do path
        a = jldopen(f -> f["timeseries/$name/$iteration"][:, :, 1], path)
        Hx = (size(a, 1) % 360) ÷ 2
        Nx = size(a, 1) - 2Hx
        Ny = 170 * (Nx ÷ 360) ÷ length(paths)
        Hy = (size(a, 2) - Ny) ÷ 2
        a[Hx .+ (1:Nx+face_x), Hy .+ (1:Ny+face_y)]
    end
    face_y || return reduce(hcat, bands)
    # Neighboring bands share the face between them: keep the southern faces of each band and the top face of the last
    return hcat([b[:, 1:end-1] for b in bands[1:end-1]]..., bands[end])
end

face_averaged_speed(u, v) = @. sqrt((u[1:end-1, :]^2 + u[2:end, :]^2) / 2 + (v[:, 1:end-1]^2 + v[:, 2:end]^2) / 2)

"""Surface ocean speed (NaN on land) and sea ice speed (NaN where ℵ < 1e-3 or on land) at full resolution."""
function speeds(directory, ocean_iteration, ice_iteration)
    uo = stitch(directory, "ocean_surface", "u", ocean_iteration; face_x = true)
    vo = stitch(directory, "ocean_surface", "v", ocean_iteration; face_y = true)
    ui = stitch(directory, "sea_ice", "u", ice_iteration; face_x = true)
    vi = stitch(directory, "sea_ice", "v", ice_iteration; face_y = true)
    ℵi = stitch(directory, "sea_ice", "ℵ", ice_iteration)

    so = face_averaged_speed(uo, vo)
    si = face_averaged_speed(ui, vi)
    land = (uo[1:end-1, :] .== 0) .& (uo[2:end, :] .== 0) .& (vo[:, 1:end-1] .== 0) .& (vo[:, 2:end] .== 0)
    so[land] .= NaN
    si[(ℵi .< minimum_ice_concentration) .| land] .= NaN
    return so, si
end

# Every `stride`-th cell, with one more column that repeats the first (the periodic seam) and one more row that holds
# the cells across the north fold, where cell i of the top row neighbors cell Nx + 1 - i.
function closed(a, stride)
    Nx, Ny = size(a)
    i, j = 1:stride:Nx, 1:stride:Ny
    b = hcat(a[i, j], a[Nx + 1 .- i, Ny])
    return vcat(b, b[1:1, :])
end

cartesian(λ, φ, r) = ((@. r * cosd(φ) * cosd(λ)), (@. r * cosd(φ) * sind(λ)), (@. r * sind(φ)))

timestamp(time) = @sprintf("%s UTC, day %.2f", Dates.format(start_date + Millisecond(round(Int, 1000time)), "yyyy-mm-dd HH:MM"), time / 86400)

"""
A figure of two globes whose colors follow the observables `so` and `si` (closed, strided speed arrays) and whose
title follows `title`.
"""
function globe_figure(Nx, Ny, stride, so, si, title)
    grid = TripolarGrid(CPU(); size = (Nx, Ny, 1), z = (0, 1), halo = (7, 7, 1))
    λ = closed(Array(λnodes(grid, Center(), Center(), Center())), stride)
    φ = closed(Array(φnodes(grid, Center(), Center(), Center())), stride)
    xo, yo, zo = cartesian(λ, φ, 1)
    xi, yi, zi = cartesian(λ, φ, 1.002) # ice just above the ocean, so the two never interleave

    # South of the grid edge: a plain gray cap with no data
    λc = range(-180, 180, length = 361)
    φc = range(-90, -80, length = 21)
    xc, yc, zc = cartesian([l for l in λc, p in φc], [p for l in λc, p in φc], 1)

    fig = Figure(size = (1200, 900))
    views = ((elevation = deg2rad(40), azimuth = deg2rad(-15)), (elevation = deg2rad(-90), azimuth = deg2rad(165)))
    plots = map(enumerate(views)) do (n, view)
        # Limits fixed to the whole sphere: with the far hemisphere dropped they would otherwise shrink unevenly
        ax = Axis3(fig[1, n]; view..., aspect = :equal, limits = (-1.01, 1.01, -1.01, 1.01, -1.01, 1.01))
        # CairoMakie paints each plot whole, in order, with no depth buffer between plots, so the far hemisphere of
        # the ice and the cap would be painted over the near hemisphere of the ocean: drop the far hemisphere.
        eye = (cos(view.elevation) * cos(view.azimuth), cos(view.elevation) * sin(view.azimuth), sin(view.elevation))
        front(x, y, z) = (far = @. x * eye[1] + y * eye[2] + z * eye[3] < 0; map(a -> ifelse.(far, NaN, a), (x, y, z)))
        ocean = surface!(ax, front(xo, yo, zo)...; color = so, colorrange = ocean_speed_range, colormap = :magma,
                         nan_color = :gray90, shading = NoShading)
        ice = surface!(ax, front(xi, yi, zi)...; color = si, colorrange = ice_speed_range, colormap = :buda10,
                       nan_color = :transparent, shading = NoShading)
        surface!(ax, front(xc, yc, zc)...; color = zeros(size(xc)), colorrange = (0, 1), colormap = [:gray55, :gray55],
                 shading = NoShading)
        hidedecorations!(ax)
        hidespines!(ax)
        (; ocean, ice)
    end

    Colorbar(fig[2, 1], plots[1].ocean, vertical = false, flipaxis = false, width = Relative(0.4),
             label = "Surface ocean speed (m s⁻¹)", labelsize = 16)
    Colorbar(fig[2, 2], plots[1].ice, vertical = false, flipaxis = false, width = Relative(0.4),
             label = "Sea ice speed (m s⁻¹), where concentration ≥ 10⁻³", labelsize = 16)
    Label(fig[0, :], title, fontsize = 18)
    colgap!(fig.layout, 1, Relative(-0.1))
    return fig
end

title_text(time) = "1/12ᵒ ocean–sea ice simulation without tides, " * timestamp(time)
