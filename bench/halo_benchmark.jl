# Wall time per step of the 1/12ᵒ coupled model on 4 GPUs with halos over NCCL or over (GPU-aware Cray) MPI:
# 30 steps at 1 minute and 20 at 5 minutes to spin up, then 200 timed steps at a fixed 10 minutes, without output.
#
#   srun -n 4 julia --project=bench bench/halo_benchmark.jl <nccl | mpi>

using MPI
using NCCL
using Oceananigans.DistributedComputations: NCCLDistributed

cd(joinpath(@__DIR__, ".."))
include("setup.jl")

communication = ARGS[1]
partition = Partition(1, 4)
arch = communication == "nccl" ? NCCLDistributed(GPU(); partition) : Distributed(GPU(); partition)
rank = MPI.Comm_rank(MPI.COMM_WORLD)

grid = global_grid(arch; cells_per_degree = 12)
model = ocean_sea_ice_model(grid; max_Δt = 10minutes, end_date = start_date + Day(5))
simulation = Simulation(model; Δt = 1minute, stop_iteration = 30)
run!(simulation)
simulation.Δt = 5minutes
simulation.stop_iteration = 50
run!(simulation)

simulation.Δt = 10minutes
step_times = Float64[]
for n in 1:200
    MPI.Barrier(MPI.COMM_WORLD)
    CUDA.synchronize()
    start = time_ns()
    time_step!(simulation)
    CUDA.synchronize()
    MPI.Barrier(MPI.COMM_WORLD)
    push!(step_times, (time_ns() - start) / 1e9)
end

rank == 0 && @info @sprintf("HALOBENCH %s: median %.3f s/step, mean %.3f s/step, min %.3f, max %.3f over %d steps at Δt = 10 min",
                            communication, median(step_times), mean(step_times), minimum(step_times), maximum(step_times),
                            length(step_times))
