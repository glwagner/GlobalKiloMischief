# Internal tides in a quiescent, horizontally uniformly stratified global ocean forced only by the astronomical
# tidal potential, after Simmons, Hallberg & Arbic (2004, Deep-Sea Research II 51, 3043–3068).
#
# The stratification is N(z) = N₀ exp(z / d) everywhere (Munk's exponential profile, close to the global mean),
# carried by temperature with a linear equation of state and constant salinity. Output is the isopycnal displacement
# ζ = b′ / N² at about 1000 m, where mode-1 internal tides are large, the vertical velocity at about 200 m, and the
# sea surface height, every hour.
#
#   julia --project internal_tides.jl <CPU | CUDA | Metal> <Float64 | Float32> <cells per degree> <Nz> <Δt in minutes> <days>

using NumericalEarth
using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: znodes
using Tidejinks
using Dates
using JLD2
using Printf
import SPICE

arch_name, FT_name = ARGS[1], ARGS[2]
cells_per_degree = parse(Float64, ARGS[3])
Nz = parse(Int, ARGS[4])
Δt = parse(Float64, ARGS[5]) * minutes
stop_days = parse(Float64, ARGS[6])

arch_name == "CUDA" && @eval using CUDA
arch_name == "Metal" && @eval using Metal
arch = arch_name == "CUDA" ? GPU() : arch_name == "Metal" ? @eval(GPU(Metal.MetalBackend())) : CPU()
FT = FT_name == "Float32" ? Float32 : Float64
Oceananigans.defaults.FloatType = FT # ocean_simulation builds its physics defaults in this type

start_date = DateTime(2015, 1, 1)
end_date = start_date + Millisecond(round(Int, stop_days * 86400e3))

cd(@__DIR__)
mkpath("output")
Tidejinks.wrangle_spice_kernels("kernels.txt", "spice")
SPICE.furnsh("kernels.txt")

Nx, Ny = round(Int, 360 * cells_per_degree), round(Int, 170 * cells_per_degree)
depth = 6000
halo = (7, 7, 4)

# The bathymetry is computed once in double precision on the CPU.
cpu_grid = TripolarGrid(CPU(); size = (Nx, Ny, 1), halo, z = (-depth, 0))
bottom_height = regrid_bathymetry(cpu_grid; minimum_depth = 20, interpolation_passes = 3, major_basins = 1)

z = ExponentialDiscretization(Nz, -depth, 0; scale = 1300, mutable = true)
underlying_grid = TripolarGrid(arch, FT; size = (Nx, Ny, Nz), halo, z)
zb = Field{Center, Center, Nothing}(underlying_grid)
set!(zb, convert(Array{FT}, interior(bottom_height)))
fill_halo_regions!(zb)
grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(zb); active_cells_map = true)

# Stratification N² = N₀² exp(2z / d), so b = N₀² d / 2 exp(2z / d), carried by T with b = g α T.
N₀, d = 5.2e-3, 1300
g, α = 9.80665, 2e-4
N²(z) = N₀^2 * exp(2z / d)
Tᵢ(x, y, z) = N₀^2 * d / 2 * exp(2z / d) / (g * α)

Φ = Tidejinks.TidalPotential(grid, (start_date, end_date); reference_frame = "IAU_EARTH")
equation_of_state = LinearEquationOfState(thermal_expansion = α, haline_contraction = 0)
# CATKE mixes the bottom boundary layer, so that bottom drag reaches the whole column on shelves.
vertical_closure = NumericalEarth.Oceans.default_ocean_closure()

# Biharmonic viscosity ν = Az² / τ removes the internal-tide energy that cascades to the grid scale. With the default
# τ = 30 days, waves two cells long decay in about 7 hours.
τ = parse(Float64, get(ENV, "BIHARMONIC_TIMESCALE", "30")) * days
@inline νhb(i, j, k, grid, τ) = Oceananigans.Operators.Azᶜᶜᶜ(i, j, k, grid)^2 / τ
ν = Field{Center, Center, Center}(grid)
set!(ν, KernelFunctionOperation{Center, Center, Center}(νhb, grid, convert(FT, τ)))
closure = (vertical_closure, HorizontalScalarBiharmonicDiffusivity(; ν))
momentum_advection = get(ENV, "MOMENTUM_ADVECTION", "default") == "none" ? nothing : WENOVectorInvariant()
free_surface = NumericalEarth.Oceans.default_free_surface(grid; fixed_Δt = Δt)

bottom_drag_coefficient = parse(Float64, get(ENV, "BOTTOM_DRAG_COEFFICIENT", "0.003"))
ocean = ocean_simulation(grid; Δt, free_surface, barotropic_potential = Φ, closure, equation_of_state, momentum_advection,
                         bottom_drag_coefficient,
                         radiative_forcing = nothing, gravitational_acceleration = g, stop_time = stop_days * days)
set!(ocean.model, T = Tᵢ, S = 35)

# Isopycnal displacement at the level nearest 1000 m, and vertical velocity at the interface nearest 200 m
zᶜ = Array(znodes(underlying_grid, Center()))
zᶠ = Array(znodes(underlying_grid, Face()))
k = argmin(abs.(zᶜ .+ 1000))
z₀ = zᶜ[k]
kʷ = argmin(abs.(zᶠ .+ 200))
T = ocean.model.tracers.T
T₀, ζ_per_T = convert(FT, Tᵢ(0, 0, z₀)), convert(FT, g * α / N²(z₀))
ζ = Field((T - T₀) * ζ_per_T; indices = (:, :, k))
w = view(ocean.model.velocities.w, :, :, kʷ)
η = ocean.model.free_surface.displacement

name = @sprintf("internal_tides_%s_%s_%gdeg_Nz%d_dt%gmin_tau%gd%s", arch_name, FT_name, 1 / cells_per_degree, Nz,
                Δt / minutes, τ / days, get(ENV, "RUN_SUFFIX", ""))
ocean.output_writers[:surface] = JLD2Writer(ocean.model, (; ζ, η, w); filename = joinpath("output", name),
                                            schedule = TimeInterval(1hours), overwrite_files = true)

wall_time = Ref(time_ns())
function progress(sim)
    elapsed = (time_ns() - wall_time[]) / 1e9
    @info @sprintf("%s: day %.2f, max|η| = %.2f m, max|ζ(%.0f m)| = %.1f m, max|u| = %.2f m/s, max|w| = %.1e m/s, wall %.1f s",
                   name, time(sim) / 86400, maximum(abs, η), z₀, maximum(abs, ζ), maximum(abs, sim.model.velocities.u),
                   maximum(abs, sim.model.velocities.w), elapsed)
    wall_time[] = time_ns()
end
add_callback!(ocean, progress, TimeInterval(1days))

run!(ocean)
