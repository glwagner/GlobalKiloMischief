# Barotropic tides in a single-level global ocean: a z⋆ shallow-water model on the 1ᵒ tripolar grid, forced only
# by the astronomical tidal potential, with quadratic bottom drag and no self-attraction and loading.
#
#   julia --project single_level.jl <architecture: CPU | CUDA | Metal> <Float64 | Float32> <Δt in minutes> <days>

using NumericalEarth
using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Tidejinks
using Dates
using JLD2
using Printf
import SPICE

arch_name, FT_name = ARGS[1], ARGS[2]
Δt = parse(Float64, ARGS[3]) * minutes
stop_days = parse(Float64, ARGS[4])
# Metal is not in this environment; add it to run on Apple GPUs.
arch_name == "CUDA" && @eval using CUDA
arch_name == "Metal" && @eval using Metal
arch = arch_name == "CUDA" ? GPU() : arch_name == "Metal" ? @eval(GPU(Metal.MetalBackend())) : CPU()
FT = FT_name == "Float32" ? Float32 : Float64
Oceananigans.defaults.FloatType = FT # ocean_simulation builds its physics defaults in this type

start_date = DateTime(2015, 1, 1)
end_date = start_date + Millisecond(round(Int, stop_days * 86400e3))

# SPICE paths are limited to 80 characters, so the kernels are listed relative to this directory.
cd(@__DIR__)
mkpath("output")
Tidejinks.wrangle_spice_kernels("kernels.txt", "spice")
SPICE.furnsh("kernels.txt")

Nx, Ny, depth = 360, 170, 6000
halo = (7, 7, 4)

# The bathymetry is computed once in double precision on the CPU.
cpu_grid = TripolarGrid(CPU(); size = (Nx, Ny, 1), halo, z = (-depth, 0))
minimum_depth = parse(Float64, get(ENV, "MINIMUM_DEPTH", "10"))
bottom_height = regrid_bathymetry(cpu_grid; minimum_depth, interpolation_passes = 3, major_basins = 1)

underlying_grid = TripolarGrid(arch, FT; size = (Nx, Ny, 1), halo, z = MutableVerticalDiscretization((-depth, 0)))
zb = Field{Center, Center, Nothing}(underlying_grid)
set!(zb, convert(Array{FT}, interior(bottom_height))) # set! does not convert Float64 arrays for Metal
fill_halo_regions!(zb)

# With one level, the default minimum fractional cell height (0.2) would deepen every column shallower than 1200 m.
grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(zb; minimum_fractional_cell_height = 1e-3))

Φ = Tidejinks.TidalPotential(grid, (start_date, end_date); reference_frame = "IAU_EARTH")
free_surface = NumericalEarth.Oceans.default_free_surface(grid; fixed_Δt = Δt)
momentum_advection = Dict("default" => WENOVectorInvariant(), "none" => nothing, "centered" => VectorInvariant())[get(ENV, "MOMENTUM_ADVECTION", "default")]
viscosity = parse(Float64, get(ENV, "HORIZONTAL_VISCOSITY", "0"))
closure = viscosity > 0 ? HorizontalScalarDiffusivity(ν = viscosity) : nothing
implicit_bottom_drag = get(ENV, "IMPLICIT_BOTTOM_DRAG", "true") == "true"

# No-slip coastal walls, keeping the immersed bottom drag that ocean_simulation would otherwise build
boundary_conditions = NamedTuple()
if get(ENV, "NO_SLIP", "false") == "true"
    O = NumericalEarth.Oceans
    drag = (μ = convert(FT, 0.003), ub = zero(FT))
    no_slip = ValueBoundaryCondition(0)
    u_bottom = O.bottom_drag_bc(O.u_immersed_drag_coefficient, O.u_immersed_bottom_drag, drag, implicit_bottom_drag)
    v_bottom = O.bottom_drag_bc(O.v_immersed_drag_coefficient, O.v_immersed_bottom_drag, drag, implicit_bottom_drag)
    boundary_conditions = (u = FieldBoundaryConditions(immersed = ImmersedBoundaryCondition(bottom = u_bottom, south = no_slip, north = no_slip)),
                           v = FieldBoundaryConditions(immersed = ImmersedBoundaryCondition(bottom = v_bottom, west = no_slip, east = no_slip)))
end
ocean = ocean_simulation(grid; Δt, free_surface, barotropic_potential = Φ, closure, momentum_advection, implicit_bottom_drag, boundary_conditions,
                         radiative_forcing = nothing, stop_time = stop_days * days)
set!(ocean.model, T = 10, S = 35)

name = @sprintf("single_level_%s_%s_dt%gmin%s", arch_name, FT_name, Δt / minutes, get(ENV, "RUN_SUFFIX", ""))
η = ocean.model.free_surface.displacement
ocean.output_writers[:η] = JLD2Writer(ocean.model, (; η); filename = joinpath("output", name),
                                      schedule = TimeInterval(1hours), array_type = Array{Float64}, overwrite_files = true)

wall_time = Ref(time_ns())
function progress(sim)
    elapsed = (time_ns() - wall_time[]) / 1e9
    @info @sprintf("%s: day %.2f, max|η| = %.3f m, max|u| = %.3f m/s, wall %.1f s",
                   name, time(sim) / 86400, maximum(abs, η), maximum(abs, sim.model.velocities.u), elapsed)
    wall_time[] = time_ns()
end
add_callback!(ocean, progress, TimeInterval(1days))

run!(ocean)
