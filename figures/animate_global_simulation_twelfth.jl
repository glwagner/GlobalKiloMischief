# Frames of the two-globe view (global_simulation_globes.jl) for every snapshot time after the start that has both
# surface ocean and sea ice output. Renders frames `first:step:last` of the time-ordered list to <frame directory>/frame_<n>.png, so that
# many processes can render disjoint ranges; slurm/animate_globes.sbatch runs them and joins the frames into an MP4.
#
#   julia --project=<env> animate_global_simulation_twelfth.jl <run directory> <frame directory> <first> <step> [stride = 2]
#   julia --project=<env> animate_global_simulation_twelfth.jl <run directory> count  # prints the number of frames

include(joinpath(@__DIR__, "global_simulation_globes.jl"))

directory = ARGS[1]
ice_snapshots = Dict(t => i for (i, t) in snapshots(directory, "sea_ice"))
# The ocean starts at rest, so the first snapshot is skipped
frames = [(i, ice_snapshots[t], t) for (i, t) in snapshots(directory, "ocean_surface") if t > 0 && haskey(ice_snapshots, t)]

if ARGS[2] == "count"
    println(length(frames))
    exit()
end

frame_directory = mkpath(ARGS[2])
first_frame, frame_step = parse(Int, ARGS[3]), parse(Int, ARGS[4])
stride = parse(Int, get(ARGS, 5, "2"))

ocean_iteration, ice_iteration, time = frames[first_frame]
so, si = speeds(directory, ocean_iteration, ice_iteration)
Nx, Ny = size(so)
so_frame = Observable(closed(so, stride))
si_frame = Observable(closed(si, stride))
title = Observable(title_text(time))
fig = globe_figure(Nx, Ny, stride, so_frame, si_frame, title)

for n in first_frame:frame_step:length(frames)
    ocean_iteration, ice_iteration, time = frames[n]
    so, si = speeds(directory, ocean_iteration, ice_iteration)
    so_frame[] = closed(so, stride)
    si_frame[] = closed(si, stride)
    title[] = title_text(time)
    save(joinpath(frame_directory, @sprintf("frame_%04d.png", n)), fig; px_per_unit = 1)
    @info "frame $n of $(length(frames)): $(timestamp(time))"
end
