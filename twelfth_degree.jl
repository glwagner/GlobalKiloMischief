# Global 1/12ᵒ ocean--sea ice simulation on 4 GPUs of one node: NCCL carries the device halo exchanges,
# MPI handles launch and scalar reductions. TripolarGrid supports y or even-x pencil partitions;
# a y partition keeps the sea ice EVP halos (2 × substeps + 3 cells) inside each 510-row subdomain.
#
#   mpiexec -n 4 julia --project twelfth_degree.jl [stop_days=30] [max_Δt_minutes=30] [cfl=1.0] [Nz=100] [surface_output_hours=6]

using MPI
using Profile  # a SIGUSR1 then prints where every thread is
using NCCL
using Oceananigans.DistributedComputations: NCCLDistributed

@time "loading packages and setup.jl" include("setup.jl")

stop_days = parse(Int, get(ARGS, 1, "30"))
max_Δt_minutes = parse(Int, get(ARGS, 2, "30"))
cfl = parse(Float64, get(ARGS, 3, "1.0"))
Nz = parse(Int, get(ARGS, 4, "100"))
surface_output_interval = parse(Int, get(ARGS, 5, "6")) * hours
max_Δt = max_Δt_minutes * minutes

arch = NCCLDistributed(GPU(); partition = Partition(1, 4))
grid = @time "building the grid" global_grid(arch; cells_per_degree = 12, Nz)
model = @time "building the model" ocean_sea_ice_model(grid; max_Δt, end_date = start_date + Day(stop_days + 1))
# INITIAL_CHECKPOINT=<directory>:<iteration> starts from a 100-level checkpoint of this grid, interpolated to Nz levels
# and relabeled to start_date, instead of from GLORYS at rest.
initial_checkpoint = get(ENV, "INITIAL_CHECKPOINT", "")
name = "twelfth_degree_dt$(max_Δt_minutes)min" * (cfl == 0.5 ? "" : "_cfl$cfl") * (Nz == 100 ? "" : "_Nz$Nz") *
       (isempty(initial_checkpoint) ? "" : "_from_checkpoint")
simulation = global_simulation(model; name, Δt = 30, max_Δt, stop_time = stop_days * days, cfl,
                               surface_output_interval)  # frames for animations

directory = joinpath(run_directory, name)
if !isempty(initial_checkpoint) && !any(contains(r"^checkpoint.*_iteration[1-9]"), readdir(directory))
    include("regrid/set_from_checkpoint.jl")
    checkpoint_directory, checkpoint_iteration = split(initial_checkpoint, ':')
    @time "interpolating $initial_checkpoint to $Nz levels" set_from_checkpoint!(model, checkpoint_directory;
        iteration = parse(Int, checkpoint_iteration), rank = MPI.Comm_rank(MPI.COMM_WORLD))
end

run!(simulation; pickup = true)
