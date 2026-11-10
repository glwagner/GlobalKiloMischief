# Dump the initial fields and three tendencies of a model built by barotropic_tides.jl or TideEye's tideeye_driver.jl,
# so the two can be compared on the CPU (analysis/diag/compare_dumps.jl). refresh_sal!(model) recomputes the SAL of the
# current η in whichever way the script does it.
using Oceananigans, JLD2, Random
using Oceananigans.Operators: Δzᶠᶜᶜ, Δzᶜᶠᶜ
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
using Oceananigans.Models.HydrostaticFreeSurfaceModels: compute_momentum_tendencies!

function dump_fields(path, model, refresh_sal!; extra...)
    grid = model.grid
    Nx, Ny = size(grid)[1:2]
    Δzu = [Δzᶠᶜᶜ(i, j, 1, grid) for i in 1:Nx, j in 1:Ny]
    Δzv = [Δzᶜᶠᶜ(i, j, 1, grid) for i in 1:Nx, j in 1:Ny]
    u, v, η = model.velocities.u, model.velocities.v, (hasproperty(model.free_surface, :displacement) ? model.free_surface.displacement : model.free_surface.η)
    T = 12.4206012 * 3600
    function tendency(t, random)
        rng = Xoshiro(1)
        for f in (u, v, η)
            set!(f, random ? 0.1 .* randn(rng, size(interior(f))...) : zeros(size(interior(f))...))
            mask_immersed_field!(f)
        end
        fill_halo_regions!((u, v, η))
        model.clock.time = t
        refresh_sal!(model)
        compute_momentum_tendencies!(model, [])
        return Array(interior(model.timestepper.Gⁿ.u, :, :, 1)), Array(interior(model.timestepper.Gⁿ.v, :, :, 1))
    end
    Gu0, Gv0 = tendency(0.0, false)
    Gu4, Gv4 = tendency(T / 4, false)
    Gur, Gvr = tendency(0.0, true)
    jldsave(path; Δzu, Δzv, Gu0, Gv0, Gu4, Gv4, Gur, Gvr, ηr = Array(interior(η, :, :, 1)),
            g = model.free_surface.gravitational_acceleration, extra...)
    @info "Dumped fields to $path"
end
