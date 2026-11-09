# A quick CPU check of tides/stratified_tides.jl without the atmosphere or sea ice: a 4ᵒ tripolar ocean with partial
# cells on a synthetic bathymetry, forced by the equilibrium potential with SAL, thickness-weighted Coriolis, the
# scale-aware wave drag (no biharmonic viscosity), run for 3 days with the harmonic fit from day 1.
#   julia --project=. tides/cpu_check.jl

using NumericalEarth
using Oceananigans
using Oceananigans.Units
using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: λnode, φnode
using Dates: Dates, DateTime
using MPI: MPI
using Printf

const start_date = DateTime(2015, 1, 1)
include(joinpath(@__DIR__, "stratified_tides.jl"))

Nx, Ny, Nz = 92, 43, 10
depth = 5000
underlying_grid = TripolarGrid(CPU(); size = (Nx, Ny, Nz), halo = (7, 7, 5), z = ExponentialDiscretization(Nz, -depth, 0; scale = 1100, mutable = true))
# Land north of 50ᵒN covers the tripolar grid's singular points (near 55ᵒN), and south of 70ᵒS
bottom_height = [ifelse(!(-70 < φnode(i, j, 1, underlying_grid, Center(), Center(), Center()) < 50), 0.0,
                        -depth * (0.3 + 0.7 * abs(sind(2λnode(i, j, 1, underlying_grid, Center(), Center(), Center()))))) for i in 1:Nx, j in 1:Ny, k in 1:1]
zb = Field{Center, Center, Nothing}(underlying_grid)
set!(zb, bottom_height)
fill_halo_regions!(zb)
grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(zb); active_cells_map = true)

max_Δt = 30minutes
ocean_kw = tidal_ocean_kw(grid, bottom_height; max_Δt, wave_drag = :scale_aware, cells_per_degree = 6)
ocean = ocean_simulation(grid; Δt = max_Δt, free_surface = NumericalEarth.Oceans.default_free_surface(grid; fixed_Δt = max_Δt),
                         radiative_forcing = nothing, stop_time = 3days, ocean_kw...)
set!(ocean.model, T = 10, S = 35)
set!(ocean.model, T = (λ, φ, z) -> 2 + 18 * exp(z / 800))  # stratified, so that N at the bottom is nonzero
set_scale_aware_wave_drag!(ocean_kw.forcing.u.rate, ocean.model, subgrid_roughness_over_2H(grid, bottom_height, 6); Δt = max_Δt)
update_self_attraction_and_loading_every_step!(ocean)
η = ocean.model.free_surface.displacement
directory = mkpath(joinpath(tempdir(), "tides_cpu_check"))
foreach(rm, filter(contains("harmonics"), readdir(directory; join = true)))
add_harmonic_fit!(ocean, η, directory; start = 1days, save_interval = 1days, first_window = 2days, minimum_span = 1days)
add_callback!(ocean, sim -> @info(@sprintf("day %.2f: max|η| %.3f m, max|u| %.3f m/s", time(sim) / day, maximum(abs, η), maximum(abs, sim.model.velocities.u))),
              TimeInterval(12hours))
run!(ocean)

fit = load(joinpath(directory, "harmonics_rank0.jld2"))
@info @sprintf("Harmonic fit: %d samples, days %.2f–%.2f, max |M2| %.3f m, max |K1| %.3f m, files %s", length(fit["times"]),
               fit["first_day"], fit["last_day"], maximum(abs, fit["Z"]["M2"]), maximum(abs, fit["Z"]["K1"]),
               join(filter(contains("harmonics"), readdir(directory)), ", "))
