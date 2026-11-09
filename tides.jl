# Stratified global ocean--sea ice simulation with astronomical tides (M2, S2, K1, O1), self-attraction and loading,
# partial cells on bathymetry regridded 11000 m deep, and an online harmonic fit of the sea surface height from day 20
# (see tides/stratified_tides.jl). One GPU, or several with MPI (a y partition). Output: surface T, S, u, v and η every
# 6 hours, hourly surface u, v and η from day 20 to day 50, sea ice daily.
#
#   julia --project tides.jl <cells per degree> <stop_days> <max_Δt_minutes> <cfl> <none | jsl | scale_aware> <χ> <biharmonic: yes | no>
#   mpiexec -n 4 julia --project tides.jl 12 ...

using MPI
using NCCL
using Oceananigans.DistributedComputations: NCCLDistributed

@time "loading packages and setup.jl" include("setup.jl")
include("tides/stratified_tides.jl")

cells_per_degree = parse(Int, ARGS[1])
stop_days = parse(Int, ARGS[2])
max_Δt_minutes = parse(Int, ARGS[3])
max_Δt = max_Δt_minutes * minutes
cfl = parse(Float64, ARGS[4])
wave_drag = Symbol(ARGS[5])
wave_drag_coefficient = parse(Float64, ARGS[6])
biharmonic = ARGS[7] == "yes"
Nz = 100

MPI.Init()
ranks = MPI.Comm_size(MPI.COMM_WORLD)
arch = ranks == 1 ? GPU() : NCCLDistributed(GPU(); partition = Partition(1, ranks))

minimum_depth = tidal_minimum_depth(Nz)
bottom_height = tidal_bottom_height(360cells_per_degree, 170cells_per_degree; minimum_depth)
depth = ceil(-minimum(bottom_height))
@info @sprintf("Tidal grid: minimum depth %.2f m, deepest cell %.0f m", minimum_depth, depth)
grid = @time "building the grid" global_grid(arch; cells_per_degree, Nz, depth, minimum_depth, bottom_height,
                                              immersed_boundary = PartialCellBottom)

ocean_kw = tidal_ocean_kw(grid, bottom_height; max_Δt, wave_drag, wave_drag_coefficient, cells_per_degree, biharmonic)
model = @time "building the model" ocean_sea_ice_model(grid; max_Δt, end_date = start_date + Day(stop_days + 1), ocean_kw)
# From the initial (GLORYS) stratification, so every segment of a run has the same drag.
wave_drag == :scale_aware && set_scale_aware_wave_drag!(ocean_kw.forcing.u.rate, model.ocean.model,
                                                        subgrid_roughness_over_2H(grid, bottom_height, cells_per_degree);
                                                        χ = wave_drag_coefficient, Δt = max_Δt)
update_self_attraction_and_loading_every_step!(model.ocean)

name = @sprintf("tides_%ddeg_dt%dmin_cfl%g_%s%g%s", cells_per_degree, max_Δt_minutes, cfl, wave_drag, wave_drag_coefficient,
                biharmonic ? "_biharmonic" : "")
simulation = global_simulation(model; name, Δt = 1minute, max_Δt, stop_time = stop_days * days, cfl,
                               surface_output_interval = 6hours, surface_fields = (:T, :S, :u, :v))

directory = joinpath(run_directory, name)
ocean = model.ocean
η = ocean.model.free_surface.displacement
hourly = SpecifiedTimes(collect(20days:1hour:50days))
ocean.output_writers[:hourly_velocities] = JLD2Writer(ocean.model, (; ocean.model.velocities.u, ocean.model.velocities.v);
                                                      dir = directory, filename = "ocean_hourly_velocities", indices = (:, :, Nz),
                                                      schedule = hourly, array_type = Array{Float32})
ocean.output_writers[:hourly_free_surface] = JLD2Writer(ocean.model, (; η); dir = directory, filename = "ocean_hourly_free_surface",
                                                        schedule = hourly, array_type = Array{Float32})
add_harmonic_fit!(simulation, η, directory)

run!(simulation; pickup = true)
