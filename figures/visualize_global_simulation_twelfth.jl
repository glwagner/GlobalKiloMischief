# Two globe views of surface ocean speed with sea ice speed on top, from the last snapshot of a distributed 1/12ᵒ run.
# Adapted from OceananigansPaper figures/visualize_global_simulation3.jl, with these changes:
#
# - Speed from face values averaged to cell centers: s = √((uᵢ² + uᵢ₊₁²)/2 + (vⱼ² + vⱼ₊₁²)/2). The original summed the
#   two faces, which made speeds √2 too large.
# - Sea ice is drawn where its concentration ℵ ≥ 1e-3; elsewhere the ocean shows through.
# - The surface mesh is closed across the periodic longitude seam and across the north fold.
# - The grid starts at 80ᵒS: south of it is a plain neutral gray cap, with no data and no label.
# - Fixed physical color ranges instead of fractions of the maxima.
#
#   julia --project=<env with CairoMakie, JLD2, Oceananigans> visualize_global_simulation_twelfth.jl <run directory> <output.png> [stride = 2]

using CairoMakie
using JLD2
using Oceananigans
using Dates
using Printf

run_directory = ARGS[1]
output_path = ARGS[2]
stride = parse(Int, get(ARGS, 3, "2"))  # plot every `stride`-th cell in each direction
start_date = DateTime(2015, 1, 1)

ocean_speed_range = (0, 1.0)   # m s⁻¹
ice_speed_range = (0, 0.4)     # m s⁻¹
minimum_ice_concentration = 1e-3

rank_files(prefix) = sort(filter(f -> startswith(f, prefix * "_rank"), readdir(run_directory)),
                          by = f -> parse(Int, match(r"rank(\d+)", f)[1]))

last_iteration(path) = jldopen(path) do file
    last(sort(filter(!=("serialized"), collect(keys(file["timeseries/t"]))), by = x -> parse(Int, x)))
end

# Stack the interiors of the ranks' bands of rows. Each array has a halo of `(size % 360) / 2` cells in x and
# `(size % 170) / 2` in y (the grid is 360c × 170c), and faces keep one more column or row than centers.
function stitch(prefix, name, iteration; face_x = false, face_y = false)
    bands = map(rank_files(prefix)) do file
        a = jldopen(joinpath(run_directory, file)) do f
            f["timeseries/$name/$iteration"][:, :, 1]
        end
        Hx, Hy = (size(a, 1) % 360) ÷ 2, (size(a, 2) % 170) ÷ 2
        Nx, Ny = size(a, 1) - 2Hx, size(a, 2) - 2Hy
        a[Hx .+ (1:Nx+face_x), Hy .+ (1:Ny+face_y)]
    end
    face_y || return reduce(hcat, bands)
    # Neighboring bands share the face between them: keep the southern faces of each band and the top face of the last
    return hcat([b[:, 1:end-1] for b in bands[1:end-1]]..., bands[end])
end

face_averaged_speed(u, v) = @. sqrt((u[1:end-1, :]^2 + u[2:end, :]^2) / 2 + (v[:, 1:end-1]^2 + v[:, 2:end]^2) / 2)

ocean_iteration = last_iteration(joinpath(run_directory, first(rank_files("ocean_surface"))))
ice_iteration = last_iteration(joinpath(run_directory, first(rank_files("sea_ice"))))
time = jldopen(f -> f["timeseries/t/$ocean_iteration"], joinpath(run_directory, first(rank_files("ocean_surface"))))
ice_time = jldopen(f -> f["timeseries/t/$ice_iteration"], joinpath(run_directory, first(rank_files("sea_ice"))))
time == ice_time || error("the last ocean snapshot ($time s) and sea ice snapshot ($ice_time s) differ")
date = start_date + Millisecond(round(Int, 1000time))
timestamp = @sprintf("%s UTC, day %.2f", Dates.format(date, "yyyy-mm-dd HH:MM"), time / 86400)
@info "Snapshot $timestamp (ocean iteration $ocean_iteration, sea ice iteration $ice_iteration)"

uo = stitch("ocean_surface", "u", ocean_iteration; face_x = true)
vo = stitch("ocean_surface", "v", ocean_iteration; face_y = true)
ui = stitch("sea_ice", "u", ice_iteration; face_x = true)
vi = stitch("sea_ice", "v", ice_iteration; face_y = true)
ℵi = stitch("sea_ice", "ℵ", ice_iteration)

so = face_averaged_speed(uo, vo)
si = face_averaged_speed(ui, vi)
Nx, Ny = size(so)
@info @sprintf("%d × %d cells; max ocean speed %.2f m/s, max ice speed %.2f m/s", Nx, Ny, maximum(so), maximum(si[ℵi .≥ minimum_ice_concentration]))

land = (uo[1:end-1, :] .== 0) .& (uo[2:end, :] .== 0) .& (vo[:, 1:end-1] .== 0) .& (vo[:, 2:end] .== 0)
so[land] .= NaN
si[(ℵi .< minimum_ice_concentration) .| land] .= NaN

grid = TripolarGrid(CPU(); size = (Nx, Ny, 1), z = (0, 1), halo = (7, 7, 1))
λ = Array(λnodes(grid, Center(), Center(), Center()))
φ = Array(φnodes(grid, Center(), Center(), Center()))

# Every `stride`-th cell, with one more column that repeats the first (the periodic seam) and one more row that holds
# the cells across the north fold, where cell i of the top row neighbors cell Nx + 1 - i.
function closed(a)
    i, j = 1:stride:Nx, 1:stride:Ny
    b = a[i, j]
    b = hcat(b, a[Nx + 1 .- i, Ny])
    return vcat(b, b[1:1, :])
end

cartesian(λ, φ, r) = (@. r * cosd(φ) * cosd(λ), @. r * cosd(φ) * sind(λ), @. r * sind(φ))

xo, yo, zo = cartesian(closed(λ), closed(φ), 1)
xi, yi, zi = cartesian(closed(λ), closed(φ), 1.002) # ice just above the ocean, so the two never interleave

# South of the grid edge: a plain gray cap with no data
λc = range(-180, 180, length = 361)
φc = range(-90, -80, length = 21)
xc, yc, zc = cartesian([l for l in λc, p in φc], [p for l in λc, p in φc], 1)

fig = Figure(size = (1200, 900))
views = ((elevation = deg2rad(40), azimuth = deg2rad(-15)), (elevation = deg2rad(-90), azimuth = deg2rad(165)))
plots = map(enumerate(views)) do (n, view)
    ax = Axis3(fig[1, n]; view..., aspect = :equal)
    ocean = surface!(ax, xo, yo, zo; color = closed(so), colorrange = ocean_speed_range, colormap = :magma,
                     nan_color = :gray90, shading = NoShading)
    ice = surface!(ax, xi, yi, zi; color = closed(si), colorrange = ice_speed_range, colormap = :buda10,
                   nan_color = :transparent, shading = NoShading)
    surface!(ax, xc, yc, zc; color = zeros(size(xc)), colorrange = (0, 1), colormap = [:gray55, :gray55], shading = NoShading)
    hidedecorations!(ax)
    hidespines!(ax)
    (; ocean, ice)
end
ocean_surface, ice_surface = plots[1].ocean, plots[1].ice

Colorbar(fig[2, 1], ocean_surface, vertical = false, flipaxis = false, width = Relative(0.4),
         label = "Surface ocean speed (m s⁻¹)", labelsize = 16)
Colorbar(fig[2, 2], ice_surface, vertical = false, flipaxis = false, width = Relative(0.4),
         label = "Sea ice speed (m s⁻¹), where concentration ≥ 10⁻³", labelsize = 16)
Label(fig[0, :], "1/12ᵒ ocean–sea ice simulation without tides, " * timestamp, fontsize = 18)
colgap!(fig.layout, 1, Relative(-0.1))

save(output_path, fig; px_per_unit = 2)
@info "Wrote $output_path"
