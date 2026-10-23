# Animate hourly sea surface height on the sphere, drawing each tripolar cell as a quadrilateral
# between its four corner (Face, Face) nodes and coloring it by the cell's η.
#
#   julia --project animate_globe.jl <run> <first day> <last day>

using Oceananigans
using Oceananigans.Grids: λnode, φnode
using CairoMakie
using CairoMakie.GeometryBasics: Point3f, GLTriangleFace
using Printf

run_name = ARGS[1]
first_day, last_day = parse.(Float64, ARGS[2:3])

η = FieldTimeSeries(joinpath(@__DIR__, "output", run_name * ".jld2"), "η"; backend = InMemory())
grid = η.grid
Nx, Ny = size(grid, 1), size(grid, 2)
wet = [!Oceananigans.ImmersedBoundaries.immersed_cell(i, j, 1, grid) for i in 1:Nx, j in 1:Ny]

function on_sphere(i, j)
    λ = deg2rad(λnode(i, j, 1, grid, Face(), Face(), Center()))
    φ = deg2rad(φnode(i, j, 1, grid, Face(), Face(), Center()))
    return Point3f(cos(φ) * cos(λ), cos(φ) * sin(λ), sin(φ))
end

# Four vertices per cell, so that each cell is colored by its own value; Makie draws each quadrilateral as two triangles
cells = [(i, j) for j in 1:Ny for i in 1:Nx]
vertices = Point3f[]
for (i, j) in cells
    append!(vertices, (on_sphere(i, j), on_sphere(i + 1, j), on_sphere(i + 1, j + 1), on_sphere(i, j + 1)))
end
faces = [face for n in eachindex(cells) for face in (GLTriangleFace(4n - 3, 4n - 2, 4n - 1), GLTriangleFace(4n - 3, 4n - 1, 4n))]

frames = findall(t -> first_day * 86400 ≤ t ≤ last_day * 86400, η.times)
cell_values(n) = [wet[i, j] ? η[n][i, j, 2] : NaN32 for (i, j) in cells]
vertex_colors(values) = repeat(values; inner = 4)

n = Observable(first(frames))
colors = @lift vertex_colors(cell_values($n))
title = @lift @sprintf("Sea surface height (m), day %.2f", η.times[$n] / 86400)

fig = Figure(size = (900, 900), backgroundcolor = :black)
Label(fig[0, 1], title; color = :white, fontsize = 22, tellwidth = false)
ax = LScene(fig[1, 1]; show_axis = false)
mesh!(ax, vertices, faces; color = colors, colormap = :balance, colorrange = (-1.5, 1.5), nan_color = :gray30, shading = NoShading)
Colorbar(fig[2, 1], colormap = :balance, colorrange = (-1.5, 1.5), vertical = false, width = Relative(0.6),
         labelcolor = :white, ticklabelcolor = :white, label = "η (m)")

cam = cam3d!(ax.scene)
distance = 3.2
function look_from!(λ, φ)
    eye = distance .* Vec3f(cosd(φ) * cosd(λ), cosd(φ) * sind(λ), sind(φ))
    update_cam!(ax.scene, eye, Vec3f(0, 0, 0), Vec3f(0, 0, 1))
end

filename = joinpath(@__DIR__, "output", run_name * "_globe.mp4")
record(fig, filename, eachindex(frames); framerate = 8) do m
    n[] = frames[m]
    look_from!(-60 + 360 * (m - 1) / length(frames), -25)
end
@info "Saved $filename"
