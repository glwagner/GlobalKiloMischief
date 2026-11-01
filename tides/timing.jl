# Where does the time go in a tidal step? Builds the 1/6ᵒ model in one configuration, takes 3 warm-up steps and times
# 20 steps of Δt = 10 min with CUDA.@sync. In configurations with spherical-harmonic SAL, also times each piece of a
# SAL update 20 times. Barotropic substeps are sized for a 60-minute cap, as in the tides smoke test.
#
#   julia --project tides/timing.jl <baseline | biharmonic | potential | sal_stage | sal_step | sal_stage_no_biharmonic>

using SpeedyTransforms, ConservativeRegridding
include(joinpath(@__DIR__, "..", "setup.jl"))
using Oceananigans.Simulations: Callback
using Statistics: median

variant = ARGS[1]
configurations = (baseline                = (tides = false, sal = false, biharmonic = false),
                  biharmonic              = (tides = false, sal = false, biharmonic = true),
                  potential               = (tides = true,  sal = false, biharmonic = true),
                  sal_stage               = (tides = true,  sal = true,  biharmonic = true),
                  sal_step                = (tides = true,  sal = true,  biharmonic = true),
                  sal_stage_no_biharmonic = (tides = true,  sal = true,  biharmonic = false))
c = configurations[Symbol(variant)]

Δt = 10minutes
grid = global_grid(GPU(); cells_per_degree = 6, minimum_depth = 20)
sal = c.sal ? SphericalHarmonicSelfAttractionAndLoading(grid; truncation = 40, nlat_half = 24) : nothing
model = ocean_sea_ice_model(grid; max_Δt = 60minutes, end_date = start_date + Day(3), tides = c.tides,
                            self_attraction_and_loading = sal, biharmonic_timescale = c.biharmonic ? 30days : nothing)
ocean = model.ocean

if variant == "sal_step"
    # Once per time step, after the step, instead of within every Runge-Kutta stage's update_state!
    pop!(ocean.callbacks, :self_attraction_and_loading)
    ocean.callbacks[:self_attraction_and_loading_per_step] =
        Callback(sim -> NumericalEarth.Oceans.compute_self_attraction_and_loading!(sim.model, sal), IterationInterval(1))
end

synchronized_seconds(f) = (CUDA.synchronize(); t = time_ns(); f(); CUDA.synchronize(); 1e-9 * (time_ns() - t))

for n in 1:3
    if n == 1 && variant == "potential"
        # Keep the PTX of every kernel compiled in the first step, for a register and spill check with ptxas -v
        kernels = mkpath(joinpath(@__DIR__, "..", "logs", "kernels-$(get(ENV, "SLURM_JOB_ID", "local"))-$variant"))
        seconds = synchronized_seconds(() -> CUDA.@device_code dir = kernels time_step!(model, Δt))
    else
        seconds = synchronized_seconds(() -> time_step!(model, Δt))
    end
    @info @sprintf("%s warm-up step %d: %.2f s", variant, n, seconds)
end
steps = [synchronized_seconds(() -> time_step!(model, Δt)) for n in 1:20]
@info @sprintf("RESULT %s: median %.3f s/step, mean %.3f, min %.3f, max %.3f over %d steps of Δt = %d min (barotropic substeps for 60 min)",
               variant, median(steps), sum(steps) / length(steps), minimum(steps), maximum(steps), length(steps), Δt / minutes)

if c.sal
    ext = Base.get_extension(NumericalEarth, :NumericalEarthSpeedyTransformsExt)
    op = sal.operator
    pᵇ = sal.bottom_pressure_anomaly
    Ψ = sal.potential
    src = op.regridder.src_temp
    pieces = (bottom_pressure = () -> NumericalEarth.Oceans.compute_bottom_pressure_anomaly!(sal, ocean.model),
              regrid = () -> begin
                  ext.copy_to_source!(src, view(interior(pᵇ), :, :, 1))
                  ext.regrid!(op.gaussian_values, op.regridder, src)
                  op.gaussian_values .= ifelse.(isfinite.(op.gaussian_values), op.gaussian_values, 0)
              end,
              transform = () -> begin
                  copyto!(op.gaussian_field.data, op.gaussian_values)
                  ext.transform!(op.coefficients, op.gaussian_field, op.transform)
                  op.coefficients.data .*= op.weights
              end,
              synthesis = () -> ext.synthesize_loading_potential!(Ψ, pᵇ.grid, op),
              halo = () -> fill_halo_regions!(Ψ),
              whole_update = () -> NumericalEarth.Oceans.compute_self_attraction_and_loading!(ocean.model, sal))
    for (name, piece) in pairs(pieces)
        piece() # compile
        times = [synchronized_seconds(piece) for n in 1:20]
        @info @sprintf("RESULT %s SAL %s: median %.4f s, max %.4f s", variant, name, median(times), maximum(times))
    end
end
