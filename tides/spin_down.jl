# Does quadratic bottom drag slow the depth-mean flow in a z⋆ ocean_simulation with a split-explicit free surface?
# In a doubly periodic box with uniform flow, du/dt = -Cᴰ |u| u / h, so u(t) = u₀ / (1 + Cᴰ u₀ t / h).
using NumericalEarth, Oceananigans, Oceananigans.Units, Printf

h, u₀, Cᴰ, Δt, stop_time = 20, 0.5, 0.003, 10minutes, 6hours
for Nz in (1, 4), implicit_bottom_drag in (true, false)
    grid = RectilinearGrid(size = (8, 8, Nz), halo = (7, 7, 4), x = (0, 100kilometers), y = (0, 100kilometers),
                           z = MutableVerticalDiscretization(range(-h, 0, length = Nz + 1)), topology = (Periodic, Periodic, Bounded))
    free_surface = SplitExplicitFreeSurface(grid; substeps = 30)
    ocean = ocean_simulation(grid; Δt, free_surface, closure = nothing, radiative_forcing = nothing, stop_time,
                             bottom_drag_coefficient = Cᴰ, implicit_bottom_drag, barotropic_potential = nothing, warn = false)
    set!(ocean.model, u = u₀, T = 10, S = 35)
    run!(ocean)
    ū = sum(interior(ocean.model.velocities.u)) / length(interior(ocean.model.velocities.u))
    @info @sprintf("Nz = %d, implicit drag %-5s: depth-mean u after 6 h = %.4f m/s (drag predicts %.4f, no drag %.4f)",
                   Nz, implicit_bottom_drag, ū, u₀ / (1 + Cᴰ * u₀ * stop_time / h), u₀)
end
