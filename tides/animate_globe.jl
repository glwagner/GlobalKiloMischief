# Animate an hourly horizontal slice on the sphere, drawing each tripolar cell as a quadrilateral between its four
# corner (Face, Face) nodes and coloring it by the cell's value. The color range is the 99th percentile of |value|
# over ocean cells in the last frame.
#
#   julia --project animate_globe.jl <run> <field: η | ζ | w> <first day> <last day>

using Oceananigans
using Oceananigans.Grids: λnode, φnode
using Oceananigans.ImmersedBoundaries: immersed_cell
using CairoMakie
using CairoMakie.GeometryBasics: Point3f, GLTriangleFace
using Printf

run_name, field_name = ARGS[1], ARGS[2]
first_day, last_day = parse.(Float64, ARGS[3:4])
labels = Dict("η" => "Sea surface height (m)", "ζ" => "Isopycnal displacement (m)", "w" => "Vertical velocity (m s⁻¹)")

c = FieldTimeSeries(joinpath(@__DIR__, "output", run_name * ".jld2"), field_name; backend = InMemory())
grid = c.grid
Nx, Ny, Nz = size(grid)
k = min(first(c.indices[3]), Nz) # the level of the slice; η lives at Nz + 1
wet = [!immersed_cell(i, j, k, grid) for i in 1:Nx, j in 1:Ny]

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

# With HIGH_PASS=true, each frame shows the departure from the centered 25-hour (about two M2 periods) mean
high_pass = get(ENV, "HIGH_PASS", "false") == "true"
window = high_pass ? 12 : 0
frames = findall(t -> first_day * 86400 ≤ t ≤ last_day * 86400, c.times)
frames = filter(n -> n - window ≥ 1 && n + window ≤ length(c.times), frames)

slice(n) = Array(interior(c[n], :, :, 1))
function cell_values(n)
    s = high_pass ? slice(n) .- sum(slice, n-window:n+window) ./ (2window + 1) : slice(n)
    return [wet[i, j] ? s[i, j] : NaN32 for (i, j) in cells]
end

ocean_values = sort(abs.(filter(!isnan, cell_values(last(frames)))))
limit = ocean_values[ceil(Int, 0.99 * length(ocean_values))]
colorrange = (-limit, limit)

n = Observable(first(frames))
colors = @lift repeat(cell_values($n); inner = 4)
title = @lift @sprintf("%s%s, day %.2f", labels[field_name], high_pass ? ", minus its 25-hour mean" : "", c.times[$n] / 86400)

fig = Figure(size = (900, 900), backgroundcolor = :black)
Label(fig[0, 1], title; color = :white, fontsize = 22, tellwidth = false)
ax = LScene(fig[1, 1]; show_axis = false)
mesh!(ax, vertices, faces; color = colors, colormap = :balance, colorrange, nan_color = :gray30, shading = NoShading)
Colorbar(fig[2, 1]; colormap = :balance, colorrange, vertical = false, width = Relative(0.6),
         labelcolor = :white, ticklabelcolor = :white, label = labels[field_name])

cam3d!(ax.scene)
function look_from!(λ, φ)
    eye = 3.2 .* Vec3f(cosd(φ) * cosd(λ), cosd(φ) * sind(λ), sind(φ))
    update_cam!(ax.scene, eye, Vec3f(0, 0, 0), Vec3f(0, 0, 1))
end

filename = joinpath(@__DIR__, "output", run_name * "_" * field_name * (high_pass ? "_high_pass" : "") * "_globe.mp4")
record(fig, filename, eachindex(frames); framerate = 8) do m
    n[] = frames[m]
    look_from!(-160 + 360 * (m - 1) / length(frames), 10)
end
@info "Saved $filename"
