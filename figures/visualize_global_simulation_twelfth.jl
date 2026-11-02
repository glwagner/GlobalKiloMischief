# Two globe views of surface ocean and sea ice speed from the last snapshot of a distributed run (see
# global_simulation_globes.jl for what is drawn and how it differs from OceananigansPaper's visualize_global_simulation3.jl).
#
#   julia --project=<env with CairoMakie, JLD2, Oceananigans> visualize_global_simulation_twelfth.jl <run directory> <output.png> [stride = 2]

include(joinpath(@__DIR__, "global_simulation_globes.jl"))

directory = ARGS[1]
output_path = ARGS[2]
stride = parse(Int, get(ARGS, 3, "2"))  # plot every `stride`-th cell in each direction

ocean_iteration, time = last(snapshots(directory, "ocean_surface"))
ice_snapshots = Dict(t => i for (i, t) in snapshots(directory, "sea_ice"))
haskey(ice_snapshots, time) || error("no sea ice snapshot at the last ocean snapshot time ($time s)")
ice_iteration = ice_snapshots[time]
@info "Snapshot $(timestamp(time)) (ocean iteration $ocean_iteration, sea ice iteration $ice_iteration)"

so, si = speeds(directory, ocean_iteration, ice_iteration)
Nx, Ny = size(so)
@info @sprintf("%d × %d cells; max ocean speed %.2f m/s, max ice speed %.2f m/s", Nx, Ny,
               maximum(filter(!isnan, so)), maximum(filter(!isnan, si)))

fig = globe_figure(Nx, Ny, stride, closed(so, stride), closed(si, stride), title_text(time))
save(output_path, fig; px_per_unit = 2)
@info "Wrote $output_path"
