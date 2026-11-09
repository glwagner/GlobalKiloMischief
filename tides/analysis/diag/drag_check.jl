# Is quadratic bottom drag in a single-layer PartialCellBottom column divided by the column depth H (correct),
# by the full 6000 m layer, or applied twice (bottom + immersed)?
using NumericalEarth, Oceananigans, Oceananigans.Units, Printf
for H in (100.0, 1000.0, 5000.0)
    ug = RectilinearGrid(CPU(); size = (16, 16, 1), x = (0, 16e5), y = (0, 16e5), z = (-6000, 0), topology = (Periodic, Periodic, Bounded), halo = (7, 7, 4))
    grid = ImmersedBoundaryGrid(ug, PartialCellBottom((x, y) -> -H); active_cells_map = true)
    Δt = 60.0
    ocean = ocean_simulation(grid; Δt, coriolis = nothing, closure = nothing, momentum_advection = nothing,
                             tracer_advection = nothing, equation_of_state = LinearEquationOfState(),
                             radiative_forcing = nothing, barotropic_potential = nothing,
                             free_surface = SplitExplicitFreeSurface(grid; substeps = 30), bottom_drag_coefficient = 0.003)
    set!(ocean.model, T = 10, S = 35, u = 1.0)
    u0 = 1.0
    time_step!(ocean.model, Δt)
    u1 = ocean.model.velocities.u[2, 2, 1]
    @printf("H = %6.0f: du/dt = %.3e, expected -c_d u²/H = %.3e, ratio %.3f; Δz = %.1f\n", H, (u1 - u0) / Δt, -0.003 / H,
            (u1 - u0) / Δt / (-0.003 / H), Oceananigans.Operators.Δzᶠᶜᶜ(2, 2, 1, grid))
end
