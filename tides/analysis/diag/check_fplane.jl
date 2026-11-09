# Free inertia-gravity wave in the single-layer split-explicit model: numerical frequency and damping vs ω² = f² + c²k².
using NumericalEarth, Oceananigans, Oceananigans.Units, CUDA, Printf, LinearAlgebra
const O = NumericalEarth.Oceans
for (f, Δtm, ts_name, sub) in ((1e-4, 10.0, :QuasiAdamsBashforth2, 0), (1e-4, 10.0, :SplitRungeKutta3, 200), (1e-4, 10.0, :SplitRungeKutta3, 30))
    Hd = 4000.0; g = 9.80665; L = 6e6; Nx = 240
    k = 2π / L * 2  # wavelength 3000 km
    grid = RectilinearGrid(GPU(); size = (Nx, 16, 1), x = (0, L), y = (0, L / 60), z = (-Hd, 0), topology = (Periodic, Periodic, Bounded), halo = (7, 7, 4))
    Δt = Δtm * 60
    fs = sub > 0 ? SplitExplicitFreeSurface(grid; substeps = sub) : SplitExplicitFreeSurface(grid; cfl = 0.7, fixed_Δt = Δt)
    ocean = ocean_simulation(grid; momentum_advection = nothing, Δt, free_surface = fs, barotropic_potential = nothing, closure = nothing,
                             bottom_drag_coefficient = 0.0, coriolis = f == 0 ? nothing : FPlane(; f), equation_of_state = LinearEquationOfState(),
                             timestepper = ts_name, radiative_forcing = nothing, gravitational_acceleration = g, stop_time = 10days, verbose = false, tracer_advection = nothing)
    set!(ocean.model, T = 10, S = 35)
    η = ocean.model.free_surface.displacement
    set!(η, (x, y, z) -> cos(k * x))
    ts = Float64[]; es = Float64[]
    x = Array(xnodes(grid, Center())); basis = cos.(k .* x)
    add_callback!(ocean, s -> (push!(ts, time(s)); push!(es, sum(Array(interior(η, :, 1, 1)) .* basis) / sum(basis .^ 2))), TimeInterval(600.0))
    run!(ocean)
    ωe = sqrt(f^2 + g * Hd * k^2)
    # Fit a e^{-dt} cos(ωt + p) + b by scanning
    best = (Inf, 0.0, 0.0, [0.0])
    for ω in ωe .* (0.9:0.0002:1.1), d in (0:0.0025:0.5) ./ 86400
        A = hcat(exp.(-d .* ts) .* cos.(ω .* ts), exp.(-d .* ts) .* sin.(ω .* ts), ones(length(ts)))
        c = A \ es; r = norm(A * c - es)
        r < best[1] && (best = (r, ω, d, c))
    end
    print("$ts_name sub=$sub "); @printf("RESULT f=%g Δt=%.1f min: ω/ωe = %.4f (ωΔt = %.3f), damping %.3f /day, Q-equivalent e-fold %.1f days, residual %.4f, per-M2-period amp loss %.2f%%\n", f, Δtm, best[2] / ωe, ωe * Δt, best[3] * 86400, 1 / (best[3] * 86400 + 1e-12), best[1] / norm(es), 100 * (1 - exp(-best[3] * 44712)))
end
