# Global 1/12ᵒ ocean--sea ice simulation on 4 GPUs of one node: NCCL carries the device halo exchanges,
# MPI handles launch and scalar reductions. TripolarGrid supports y or even-x pencil partitions;
# a y partition keeps the sea ice EVP halos (2 × substeps + 3 cells) inside each 510-row subdomain.
#
#   mpiexec -n 4 julia --project twelfth_degree.jl [stop_days=30] [max_Δt_minutes=10] [cfl=0.5]

using MPI
using Profile  # a SIGUSR1 then prints where every thread is
using NCCL
using Oceananigans.DistributedComputations: NCCLDistributed

@time "loading packages and setup.jl" include("setup.jl")

stop_days = parse(Int, get(ARGS, 1, "30"))
max_Δt_minutes = parse(Int, get(ARGS, 2, "10"))
cfl = parse(Float64, get(ARGS, 3, "0.5"))
max_Δt = max_Δt_minutes * minutes

arch = NCCLDistributed(GPU(); partition = Partition(1, 4))
grid = @time "building the grid" global_grid(arch; cells_per_degree = 12)
model = @time "building the model" ocean_sea_ice_model(grid; max_Δt, end_date = start_date + Day(stop_days + 1))
simulation = global_simulation(model; name = "twelfth_degree_dt$(max_Δt_minutes)min" * (cfl == 0.5 ? "" : "_cfl$cfl"), Δt = 30, max_Δt, stop_time = stop_days * days, cfl,
                               surface_output_interval = 6hours)  # frames for animations

run!(simulation; pickup = true)
