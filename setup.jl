# Shared configuration for the global ocean--sea ice simulations on a TripolarGrid:
# GLORYS initial conditions, ETOPO bathymetry, ERA5 atmosphere and radiation, JRA55 runoff.

using NumericalEarth
using Oceananigans
using Oceananigans.Units
using Oceananigans.Advection: cell_advection_timescale
using Oceananigans.DistributedComputations: @root
using Oceananigans.Fields: ZeroField
using ClimaSeaIce.Rheologies: ElastoViscoPlasticRheology
using ClimaSeaIce.SeaIceDynamics: SplitExplicitSolver
using CUDA
using Dates
using Printf
using Statistics

const start_date = DateTime(2015, 1, 1)
const run_directory = get(ENV, "RUN_DIRECTORY", ".")

"""
    global_grid(arch; cells_per_degree, Nz = 100, depth = 6000)

Tripolar grid from 80°S to the north pole at `1/cells_per_degree` degree resolution, with `Nz` z⋆ levels
(1.4 m at the surface, 320 m at 6000 m depth) and ETOPO2022 bathymetry.
"""
function global_grid(arch; cells_per_degree, Nz = 100, depth = 6000, halo = (5, 5, 4))
    Nx = 360 * cells_per_degree
    Ny = 170 * cells_per_degree
    z = ExponentialDiscretization(Nz, -depth, 0; scale = 1100, mutable = true)
    underlying_grid = TripolarGrid(arch; size = (Nx, Ny, Nz), halo, z)

    # Two basins keeps the Mediterranean even if the Strait of Gibraltar is closed at this resolution.
    bottom_height = regrid_bathymetry(underlying_grid; minimum_depth = 10, interpolation_passes = 3, major_basins = 2)

    return ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height); active_cells_map = true)
end

"""
    ocean_sea_ice_model(grid; max_Δt, end_date, sea_ice_substeps = 120)

ERA5-forced ocean--sea ice model initialized from GLORYS at `start_date` and from rest.
The ocean uses the `ocean_simulation` defaults: `WENOVectorInvariant` momentum advection, `WENO(order=7)`
tracer advection (both vertically implicit where the vertical CFL exceeds 0.5) and CATKE vertical mixing only.
The barotropic substep count is sized for `max_Δt`, so any `Δt ≤ max_Δt` is stable.
"""
function ocean_sea_ice_model(grid; max_Δt, end_date, sea_ice_substeps = 120)
    arch = architecture(grid)

    land = JRA55PrescribedLand(grid; dataset = MultiYearJRA55(), start_date, end_date)

    free_surface = NumericalEarth.Oceans.default_free_surface(grid; fixed_Δt = max_Δt)
    ocean = ocean_simulation(grid; free_surface, river_routing = land.river_routing)

    # mEVP converges only if the substep count is at least the relaxation parameter.
    rheology = ElastoViscoPlasticRheology(max_relaxation_parameter = sea_ice_substeps)
    solver = SplitExplicitSolver(grid; substeps = sea_ice_substeps)
    dynamics = NumericalEarth.SeaIces.sea_ice_dynamics(grid, ocean; rheology, solver)
    sea_ice = sea_ice_simulation(grid, ocean; advection = WENO(order = 7), dynamics)

    glorys = MetadataSet(:temperature, :salinity, :sea_ice_thickness, :sea_ice_concentration;
                         dataset = GLORYSDaily(), date = start_date)
    set!(ocean.model, glorys)
    set!(sea_ice.model, glorys)

    atmosphere = ERA5PrescribedAtmosphere(arch; start_date, end_date, time_indices_in_memory = 48)
    ocean_surface = SurfaceRadiationProperties(albedo = LatitudeDependentAlbedo())
    radiation = ERA5PrescribedRadiation(arch; start_date, end_date, ocean_surface, time_indices_in_memory = 48)

    return OceanSeaIceModel(ocean, sea_ice; atmosphere, land, radiation)
end

# Vertical advection is implicit where it would limit Δt, so only horizontal advection sets the time step.
horizontal_advection_timescale(model) =
    cell_advection_timescale(model.ocean.model.grid, (; model.ocean.model.velocities.u,
                                                        model.ocean.model.velocities.v,
                                                        w = ZeroField()))

"""
    global_simulation(model; name, Δt, max_Δt, stop_time)

Coupled simulation with an adaptive time step, a progress message every 50 iterations, daily surface output,
and a checkpoint every 10 days so that a later run with a longer `stop_time` picks up where this one ended.
"""
function global_simulation(model; name, Δt, max_Δt, stop_time, checkpoint_interval = 10days)
    simulation = Simulation(model; Δt, stop_time)
    directory = mkpath(joinpath(run_directory, name))

    wizard = TimeStepWizard(; cfl = 0.5, max_Δt, max_change = 1.05, cell_advection_timescale = horizontal_advection_timescale)
    add_callback!(simulation, wizard, IterationInterval(10))
    add_callback!(simulation, progress, IterationInterval(50))

    ocean = model.ocean
    sea_ice = model.sea_ice
    Nz = size(ocean.model.grid, 3)
    surface_outputs = merge(ocean.model.tracers, ocean.model.velocities)
    ocean.output_writers[:surface] = JLD2Writer(ocean.model, surface_outputs; dir = directory,
                                                filename = "ocean_surface", indices = (:, :, Nz),
                                                schedule = TimeInterval(1day), array_type = Array{Float32},
                                                including = [:grid])

    ocean.output_writers[:free_surface] = JLD2Writer(ocean.model, (; η = ocean.model.free_surface.displacement);
                                                     dir = directory, filename = "ocean_free_surface",
                                                     schedule = TimeInterval(1day), array_type = Array{Float32})

    sea_ice_outputs = merge((; h = sea_ice.model.ice_thickness, ℵ = sea_ice.model.ice_concentration),
                            sea_ice.model.velocities)
    sea_ice.output_writers[:surface] = JLD2Writer(sea_ice.model, sea_ice_outputs; dir = directory,
                                                  filename = "sea_ice", schedule = TimeInterval(1day),
                                                  array_type = Array{Float32})

    simulation.output_writers[:checkpointer] = Checkpointer(model; dir = directory, cleanup = true,
                                                            schedule = TimeInterval(checkpoint_interval))

    return simulation
end

const wall_clock = Ref(time_ns())
const last_progress = Ref((iteration = 0, time = 0.0))

function progress(sim)
    ocean = sim.model.ocean.model
    sea_ice = sim.model.sea_ice.model
    u, v, w = ocean.velocities
    T = ocean.tracers.T
    e = ocean.tracers.e

    elapsed = 1e-9 * (time_ns() - wall_clock[])
    steps = iteration(sim) - last_progress[].iteration
    simulated = time(sim) - last_progress[].time
    sypd = simulated / elapsed / 365  # simulated days per wall-clock second / 365 = simulated years per day

    msg = @sprintf("%s (iter %d), Δt: %s, wall/step: %s, SYPD: %.2f", Date(start_date + Millisecond(round(Int, 1000 * time(sim)))),
                   iteration(sim), prettytime(sim.Δt), prettytime(elapsed / max(steps, 1)), sypd)
    msg *= @sprintf(", max|u|: (%.2f, %.2f, %.1e) m s⁻¹", maximum(abs, u), maximum(abs, v), maximum(abs, w))
    msg *= @sprintf(", extrema(T): (%.2f, %.2f) ᵒC, max(e): %.1e m² s⁻², max(hᵢ): %.2f m",
                    minimum(T), maximum(T), maximum(e), maximum(sea_ice.ice_thickness))
    msg *= @sprintf(", GPU memory: %.1f GiB", (CUDA.total_memory() - CUDA.available_memory()) / 2^30)

    @root @info msg

    wall_clock[] = time_ns()
    last_progress[] = (iteration = iteration(sim), time = time(sim))

    return nothing
end
