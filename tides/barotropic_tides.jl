# Barotropic tides in a single-layer global ocean of uniform density forced by the astronomical tidal potential,
# for comparing treatments of self-attraction and loading (SAL) with each other and with altimetric atlases.
# The grid, bathymetry, bottom drag and biharmonic viscosity are as in internal_tides.jl. There is no wave drag,
# so the tides are expected to be too large. Output is the sea surface height every hour.
#
#   SELF_ATTRACTION_AND_LOADING=<none | scalar | spherical_harmonic | convolution> \
#       julia --project barotropic_tides.jl <CPU | CUDA> <cells per degree> <Δt in minutes> <days>
#
# The convolution (Chen, Wang, Arbic & Krasny 2026) is summed directly in CONVOLUTION_PRECISION (default Float32).

using NumericalEarth
using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Tidejinks
using Dates
using JLD2
using Printf
import SPICE

arch_name = ARGS[1]
cells_per_degree = parse(Float64, ARGS[2])
Δt = parse(Float64, ARGS[3]) * minutes
stop_days = parse(Float64, ARGS[4])

arch_name == "CUDA" && @eval using CUDA
arch = arch_name == "CUDA" ? GPU() : CPU()

start_date = DateTime(2015, 1, 1)
end_date = start_date + Millisecond(round(Int, stop_days * 86400e3))

cd(@__DIR__)
mkpath("output")
Tidejinks.wrangle_spice_kernels("kernels.txt", "spice")
SPICE.furnsh("kernels.txt")

Nx, Ny = round(Int, 360 * cells_per_degree), round(Int, 170 * cells_per_degree)
depth = 6000
halo = (7, 7, 4)

cpu_grid = TripolarGrid(CPU(); size = (Nx, Ny, 1), halo, z = (-depth, 0))
bottom_height = regrid_bathymetry(cpu_grid; minimum_depth = 20, interpolation_passes = 3, major_basins = 1)

underlying_grid = TripolarGrid(arch; size = (Nx, Ny, 1), halo, z = (-depth, 0))
zb = Field{Center, Center, Nothing}(underlying_grid)
set!(zb, Array(interior(bottom_height)))
fill_halo_regions!(zb)
grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(zb); active_cells_map = true)

Φ = Tidejinks.TidalPotential(grid, (start_date, end_date); reference_frame = "IAU_EARTH")
g = 9.80665

τ = parse(Float64, get(ENV, "BIHARMONIC_TIMESCALE", "30")) * days
@inline νhb(i, j, k, grid, τ) = Oceananigans.Operators.Azᶜᶜᶜ(i, j, k, grid)^2 / τ
ν = Field{Center, Center, Center}(grid)
set!(ν, KernelFunctionOperation{Center, Center, Center}(νhb, grid, τ))
closure = HorizontalScalarBiharmonicDiffusivity(; ν)
free_surface = NumericalEarth.Oceans.default_free_surface(grid; fixed_Δt = Δt)

loading = get(ENV, "SELF_ATTRACTION_AND_LOADING", "none")
self_attraction_and_loading = if loading == "scalar"
    ScalarSelfAttractionAndLoading(parse(Float64, get(ENV, "SCALAR_LOADING_COEFFICIENT", "0.1")))
elseif loading == "spherical_harmonic"
    @eval using SpeedyTransforms, ConservativeRegridding
    truncation = parse(Int, get(ENV, "LOADING_TRUNCATION", "40"))
    SphericalHarmonicSelfAttractionAndLoading(grid; truncation, nlat_half = max(24, cld(truncation + 1, 2) + 1))
elseif loading == "convolution"
    FT = get(ENV, "CONVOLUTION_PRECISION", "Float32") == "Float64" ? Float64 : Float32
    ConvolutionSelfAttractionAndLoading(grid; FT)
else
    nothing
end
@info "Self-attraction and loading: $(isnothing(self_attraction_and_loading) ? "none" : summary(self_attraction_and_loading))"

bottom_drag_coefficient = parse(Float64, get(ENV, "BOTTOM_DRAG_COEFFICIENT", "0.003"))
ocean = ocean_simulation(grid; Δt, free_surface, barotropic_potential = Φ, closure, bottom_drag_coefficient,
                         self_attraction_and_loading, equation_of_state = LinearEquationOfState(),
                         radiative_forcing = nothing, gravitational_acceleration = g, stop_time = stop_days * days)
set!(ocean.model, T = 10, S = 35)

η = ocean.model.free_surface.displacement
name = @sprintf("barotropic_tides_%gdeg_dt%gmin_%s", 1 / cells_per_degree, Δt / minutes, loading)
ocean.output_writers[:surface] = JLD2Writer(ocean.model, (; η); filename = joinpath("output", name),
                                            schedule = TimeInterval(1hours), overwrite_files = true)

wall_time = Ref(time_ns())
function progress(sim)
    elapsed = (time_ns() - wall_time[]) / 1e9
    @info @sprintf("%s: day %.2f, max|η| = %.2f m, max|u| = %.2f m/s, wall %.1f s",
                   name, time(sim) / 86400, maximum(abs, η), maximum(abs, sim.model.velocities.u), elapsed)
    flush(stderr)
    wall_time[] = time_ns()
end
add_callback!(ocean, progress, TimeInterval(1days))

run!(ocean)
