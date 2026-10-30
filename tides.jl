# Global ocean--sea ice simulation with astronomical tides and self-attraction and loading, for movies of eddies and
# tides together: hourly surface velocity and sea surface height. One GPU, or several with MPI (a y partition).
# The minimum depth is 20 m (10 m elsewhere), so that the tidal range cannot drain the shallowest columns.
#
#   julia --project tides.jl <cells per degree> <stop_days> <max_Δt_minutes> <cfl> <none | scalar | harmonic>
#   mpiexec -n 4 julia --project tides.jl 12 ...

using MPI
using NCCL
using Oceananigans.DistributedComputations: NCCLDistributed
using SpeedyTransforms, ConservativeRegridding  # load SphericalHarmonicSelfAttractionAndLoading

@time "loading packages and setup.jl" include("setup.jl")

cells_per_degree = parse(Int, ARGS[1])
stop_days = parse(Int, ARGS[2])
max_Δt_minutes = parse(Int, ARGS[3])
max_Δt = max_Δt_minutes * minutes
cfl = parse(Float64, ARGS[4])
loading = ARGS[5]

MPI.Init()
ranks = MPI.Comm_size(MPI.COMM_WORLD)
arch = ranks == 1 ? GPU() : NCCLDistributed(GPU(); partition = Partition(1, ranks))
grid = @time "building the grid" global_grid(arch; cells_per_degree, minimum_depth = 20)

self_attraction_and_loading = if loading == "scalar"
    ScalarSelfAttractionAndLoading(0.1)
elseif loading == "harmonic"
    SphericalHarmonicSelfAttractionAndLoading(grid; truncation = 40, nlat_half = 24)  # as in LLC4320v2
else
    nothing
end

model = @time "building the model" ocean_sea_ice_model(grid; max_Δt, end_date = start_date + Day(stop_days + 1), tides = true,
                                                      self_attraction_and_loading, biharmonic_timescale = 30days)
name = @sprintf("tides_%ddeg_dt%dmin_cfl%g_%s", cells_per_degree, max_Δt_minutes, cfl, loading)
simulation = global_simulation(model; name, Δt = 1minute, max_Δt, stop_time = stop_days * days, cfl,
                               surface_output_interval = 1hour, surface_fields = (:u, :v))

run!(simulation; pickup = true)
