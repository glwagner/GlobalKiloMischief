# Shared configuration for the global ocean--sea ice simulations on a TripolarGrid:
# GLORYS initial conditions, ETOPO bathymetry, ERA5 atmosphere and radiation, JRA55 runoff.

using NumericalEarth
using Oceananigans
using Oceananigans.Units
using Oceananigans.Architectures: architecture
using Oceananigans.Advection: cell_advection_timescale
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.DistributedComputations: @root
using Oceananigans.Fields: ZeroField
using Oceananigans.Grids: λnode, φnode
using ClimaSeaIce: IncrementalRemapping
using ClimaSeaIce.Rheologies: ElastoViscoPlasticRheology
using ClimaSeaIce.SeaIceDynamics: SplitExplicitSolver
using CUDA
using CDSAPI           # loads the ERA5 download extension
using CopernicusMarine  # loads the GLORYS download extension
using Dates: Dates, Date, DateTime, Day, Millisecond
using Printf
using Statistics

const start_date = DateTime(2015, 1, 1)
const run_directory = get(ENV, "RUN_DIRECTORY", ".")

@inline function squared_angular_distance(i, j, grid, λ₀, φ₀)
    λ = λnode(i, j, 1, grid, Center(), Center(), Center())
    φ = φnode(i, j, 1, grid, Center(), Center(), Center())
    Δλ = mod(λ - λ₀ + 180, 360) - 180
    return (Δλ * cosd(φ₀))^2 + (φ - φ₀)^2
end

nearest_cell(grid, (λ₀, φ₀)) =
    argmin(((i, j),) -> squared_angular_distance(i, j, grid, λ₀, φ₀), Iterators.product(1:size(grid, 1), 1:size(grid, 2)))

"""
    connect_basins!(bottom_height, from, to; depth)

Ensure that the ocean cells nearest `from = (λ, φ)` and `to` belong to the same basin, meaning the same region of
wet cells connected through shared faces (diagonal contact carries no flow on a C grid). If they belong to different
basins, the cells along the path between them that needs the least total deepening are deepened to `depth`. The path
is searched within the index box spanned by the two cells, padded by its own size, so it stays local to the passage
and the result does not depend on resolution. Returns the number of deepened cells.
"""
function connect_basins!(bottom_height, from, to; depth)
    grid = bottom_height.grid
    z = view(interior(bottom_height), :, :, 1)

    i₁, j₁ = nearest_cell(grid, from)
    i₂, j₂ = nearest_cell(grid, to)

    basins = NumericalEarth.Bathymetry.ImageMorphology.label_components(z .< 0)
    basins[i₁, j₁] == basins[i₂, j₂] != 0 && return 0

    pad = max(2, abs(i₂ - i₁), abs(j₂ - j₁))
    is = max(1, min(i₁, i₂) - pad):min(size(grid, 1), max(i₁, i₂) + pad)
    js = max(1, min(j₁, j₂) - pad):min(size(grid, 2), max(j₁, j₂) + pad)
    window = view(z, is, js)

    # Least-deepening path by Bellman-Ford relaxation: the cost of entering a cell is the deepening it needs.
    deepening = @. max(0, window + depth)
    start = CartesianIndex(i₁ - first(is) + 1, j₁ - first(js) + 1)
    stop  = CartesianIndex(i₂ - first(is) + 1, j₂ - first(js) + 1)
    cost = fill(Inf, size(window))
    previous = fill(start, size(window))
    cost[start] = deepening[start]
    faces = (CartesianIndex(1, 0), CartesianIndex(-1, 0), CartesianIndex(0, 1), CartesianIndex(0, -1))

    relaxed = true
    while relaxed
        relaxed = false
        for c in CartesianIndices(cost), d in faces
            n = c + d
            checkbounds(Bool, cost, n) || continue
            if cost[c] + deepening[n] < cost[n]
                cost[n] = cost[c] + deepening[n]
                previous[n] = c
                relaxed = true
            end
        end
    end

    deepened = 0
    c = stop
    while true
        if deepening[c] > 0
            window[c] = -depth
            deepened += 1
        end
        c == start && break
        c = previous[c]
    end

    return deepened
end

"""
    global_bottom_height(Nx, Ny; halo, gibraltar_sill_depth = 280)

ETOPO2022 bottom height on the whole `Nx × Ny` tripolar grid, computed on the CPU, keeping only the world ocean.
Before minor basins are removed, the Mediterranean is connected to the Atlantic through the Strait of Gibraltar if
the regridded bathymetry has closed it (at 1/6ᵒ a single land cell does), carving a channel `gibraltar_sill_depth`
deep (the Camarinal Sill). Every rank of a distributed run builds the same field, so the regridding must already be
cached (`download_data.jl` does this).
"""
function global_bottom_height(Nx, Ny; halo, gibraltar_sill_depth = 280)
    grid = TripolarGrid(CPU(); size = (Nx, Ny, 1), halo, z = (-6000, 0))
    bottom_height = regrid_bathymetry(grid; minimum_depth = 10, interpolation_passes = 3, major_basins = Inf)

    gulf_of_cadiz = (-6.8, 35.9)
    alboran_sea = (-4.6, 36.0)
    deepened = connect_basins!(bottom_height, gulf_of_cadiz, alboran_sea; depth = gibraltar_sill_depth)
    deepened > 0 && @info "Opened the Strait of Gibraltar on the $Nx × $Ny grid by deepening $deepened cells"

    NumericalEarth.Bathymetry.remove_minor_basins!(bottom_height, 1)

    return bottom_height
end

"""
    global_grid(arch; cells_per_degree, Nz = 100, depth = 6000)

Tripolar grid from 80°S to the north pole at `1/cells_per_degree` degree resolution, with `Nz` z⋆ levels
(1.4 m at the surface, 320 m at 6000 m depth) and ETOPO2022 bathymetry (see [`global_bottom_height`](@ref)).
"""
function global_grid(arch; cells_per_degree, Nz = 100, depth = 6000, halo = (7, 7, 4))
    Nx = 360 * cells_per_degree
    Ny = 170 * cells_per_degree
    z = ExponentialDiscretization(Nz, -depth, 0; scale = 1100, mutable = true)
    underlying_grid = TripolarGrid(arch; size = (Nx, Ny, Nz), halo, z)

    # `set!` partitions a global host array across the ranks of a distributed grid.
    bottom_height = Field{Center, Center, Nothing}(underlying_grid)
    set!(bottom_height, Array(interior(global_bottom_height(Nx, Ny; halo))))
    fill_halo_regions!(bottom_height)

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
    # Incremental remapping transports ice volume and concentration together and requires forward Euler.
    # WENO reconstructs them separately, which piles thick ice into coastal cells while ℵ < 1.
    sea_ice = sea_ice_simulation(grid, ocean; advection = IncrementalRemapping(), timestepper = :ForwardEuler, dynamics)

    glorys = MetadataSet(:temperature, :salinity, :sea_ice_thickness, :sea_ice_concentration;
                         dataset = GLORYSDaily(), date = start_date)
    set!(ocean.model, glorys)
    set!(sea_ice.model, glorys)

    # GLORYS marks ice-free ocean as missing rather than zero ice.
    # TODO: move into NumericalEarth's GLORYS sea ice loading.
    for field in (sea_ice.model.ice_thickness, sea_ice.model.ice_concentration)
        parent(field) .= ifelse.(isnan.(parent(field)), 0, parent(field))
    end

    atmosphere = ERA5PrescribedAtmosphere(arch; start_date, end_date, time_indices_in_memory = 48,
                                          pressure = :mean_sea_level_pressure)
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
    msg *= @sprintf(", extrema(T): (%.2f, %.2f) ᵒC, max(e): %.1e m² s⁻², max(hᵢℵ): %.2f m",
                    minimum(T), maximum(T), maximum(e), maximum(sea_ice.ice_thickness * sea_ice.ice_concentration))
    CUDA.functional() && (msg *= @sprintf(", GPU memory: %.1f GiB", (CUDA.total_memory() - CUDA.free_memory()) / 2^30))

    @root @info msg
    flush(stderr)  # batch output is block-buffered, and a crash would discard it

    wall_clock[] = time_ns()
    last_progress[] = (iteration = iteration(sim), time = time(sim))

    return nothing
end
