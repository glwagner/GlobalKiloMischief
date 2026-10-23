# Animate surface current speed from the daily `ocean_surface.jld2` output of one or more run directories, in order.
# Snapshots are deduplicated by time, so segments that redid days after a restart contribute each day once.
#
#   julia --project=<env with CairoMakie and JLD2> animate_surface_speed.jl <output.mp4> <run directory>...

using CairoMakie
using JLD2
using Dates
using Printf

output_path = ARGS[1]
run_directories = ARGS[2:end]
start_date = DateTime(2015, 1, 1)

# (time, file, iteration) of every snapshot, latest file winning for a repeated time
snapshots = Dict{Float64, Tuple{String, String}}()
for directory in run_directories
    path = joinpath(directory, "ocean_surface.jld2")
    jldopen(path) do file
        for iteration in keys(file["timeseries/t"])
            iteration == "serialized" && continue
            snapshots[file["timeseries/t/$iteration"]] = (path, iteration)
        end
    end
end
times = sort(collect(keys(snapshots)))
@info "$(length(times)) snapshots from $(Date(start_date + Second(round(Int, first(times))))) to $(Date(start_date + Second(round(Int, last(times)))))"

# Speed at cell centers from face velocities, NaN on land (where every face velocity is zero). The grids have
# 360c × 170c cells, so whatever the size exceeds a multiple of 360 or 170 by is the halo on both sides.
function surface_speed(path, iteration)
    u, v = jldopen(path) do file
        (file["timeseries/u/$iteration"][:, :, 1], file["timeseries/v/$iteration"][:, :, 1])
    end
    Hx, Hy = (size(u, 1) % 360) ÷ 2, (size(u, 2) % 170) ÷ 2
    Nx, Ny = size(u, 1) - 2Hx, size(u, 2) - 2Hy
    i, j = Hx .+ (1:Nx), Hy .+ (1:Ny)
    east = Hx > 0 ? i .+ 1 : vcat(2:Nx, 1)
    north = min.(j .+ 1, size(v, 2))
    ū = (u[i, j] .+ u[east, j]) ./ 2
    v̄ = (v[i, j] .+ v[i, north]) ./ 2
    s = sqrt.(ū .^ 2 .+ v̄ .^ 2)
    return ifelse.(s .== 0, NaN32, s)
end

n = Observable(1)
speed = @lift surface_speed(snapshots[times[$n]]...)
title = @lift @sprintf("Surface speed, %s", Date(start_date + Second(round(Int, times[$n]))))

Nx, Ny = size(speed[])
fig = Figure(size = (1600, round(Int, 1600 * Ny / Nx) + 80), backgroundcolor = :black)
ax = Axis(fig[1, 1]; title, titlecolor = :white, aspect = DataAspect(), backgroundcolor = :gray20)
hidedecorations!(ax)
hidespines!(ax)
hm = heatmap!(ax, speed; colormap = :deep, colorrange = (0, 1.2), nan_color = :gray20)
Colorbar(fig[1, 2], hm; label = "m s⁻¹", labelcolor = :white, tickcolor = :white, ticklabelcolor = :white)

save(replace(output_path, ".mp4" => "_first.png"), fig)
record(fig, output_path, eachindex(times); framerate = 12) do i
    n[] = i
end
save(replace(output_path, ".mp4" => "_last.png"), fig)
@info "Wrote $output_path"
