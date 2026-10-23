# Harmonic analysis of hourly sea surface height from single_level.jl, comparing runs against a reference run.
#
#   julia --project harmonic_analysis.jl <reference run> <other runs...>
#
# where each run is an output file name without its extension, e.g. single_level_CPU_Float64_dt10min.

using Oceananigans
using Oceananigans.Grids: λnodes, φnodes
using Oceananigans.Operators: Azᶜᶜᶜ
using CairoMakie
using LinearAlgebra
using Printf

# Periods in hours
constituents = (M2 = 12.4206012, S2 = 12.0, K1 = 23.9344696, O1 = 25.8193417)
analysis_window = (10, 30) # days

"""
Least-squares fit of η(t) = c + Σₖ (aₖ cos ωₖt + bₖ sin ωₖt) at every cell, returning
the complex amplitudes Zₖ = aₖ - i bₖ so that η ≈ Re(Zₖ exp(iωₖt)).
"""
function harmonic_fit(filename)
    η = FieldTimeSeries(joinpath(@__DIR__, "output", filename * ".jld2"), "η"; backend = InMemory())
    t = η.times
    window = findall(τ -> analysis_window[1] * 86400 ≤ τ ≤ analysis_window[2] * 86400, t)
    ω = [2π / (3600 * T) for T in values(constituents)]
    τ = t[window]
    A = hcat(ones(length(τ)), [cos(ωₖ * τⁿ) for τⁿ in τ, ωₖ in ω], [sin(ωₖ * τⁿ) for τⁿ in τ, ωₖ in ω])
    Nx, Ny = size(η.grid, 1), size(η.grid, 2)
    data = reshape(permutedims(Array(interior(η, :, :, 1, window)), (3, 1, 2)), length(window), :)
    coefficients = A \ data
    K = length(ω)
    Z = [reshape(coefficients[1 + k, :] .- im .* coefficients[1 + K + k, :], Nx, Ny) for k in 1:K]
    return NamedTuple{keys(constituents)}(Z), η.grid
end

reference_name, other_names = ARGS[1], ARGS[2:end]
reference, grid = harmonic_fit(reference_name)
Nx, Ny = size(grid, 1), size(grid, 2)
area = [Azᶜᶜᶜ(i, j, 1, grid) for i in 1:Nx, j in 1:Ny]
wet = abs.(reference.M2) .> 0
rms(x) = sqrt(sum(area[wet] .* abs2.(x[wet])) / sum(area[wet]))

@printf("%s: area-weighted rms amplitude M2 %.3f m, S2 %.3f m, K1 %.3f m, O1 %.3f m; max M2 %.2f m\n",
        reference_name, rms(reference.M2), rms(reference.S2), rms(reference.K1), rms(reference.O1), maximum(abs, reference.M2))

for name in other_names
    other, _ = harmonic_fit(name)
    for c in keys(constituents)
        @printf("%s vs reference, %s: rms |ΔZ| = %.2e m (%.1e of rms amplitude), max |ΔZ| = %.2e m\n",
                name, c, rms(other[c] .- reference[c]), rms(other[c] .- reference[c]) / rms(reference[c]),
                maximum(abs, (other[c] .- reference[c])[wet]))
    end
end

# M2 amplitude and phase of the reference run
λ = Array(λnodes(grid, Center(), Center(), Center()))[1:Nx, 1:Ny]
φ = Array(φnodes(grid, Center(), Center(), Center()))[1:Nx, 1:Ny]
amplitude = ifelse.(wet, abs.(reference.M2), NaN)
phase = ifelse.(wet, rad2deg.(angle.(reference.M2)), NaN)

fig = Figure(size = (1100, 900))
ax = Axis(fig[1, 1]; title = "M2 amplitude (m), days $(analysis_window[1])–$(analysis_window[2]): $reference_name", xlabel = "Longitude (ᵒ)", ylabel = "Latitude (ᵒ)")
sf = surface!(ax, λ, φ, zeros(Nx, Ny); color = amplitude, colormap = :viridis, colorrange = (0, 1.5), shading = NoShading)
Colorbar(fig[1, 2], sf)
ax = Axis(fig[2, 1]; title = "M2 phase (ᵒ, relative to the start of the run)", xlabel = "Longitude (ᵒ)", ylabel = "Latitude (ᵒ)")
sf = surface!(ax, λ, φ, zeros(Nx, Ny); color = phase, colormap = :romaO, colorrange = (-180, 180), shading = NoShading)
Colorbar(fig[2, 2], sf)
save(joinpath(@__DIR__, "output", reference_name * "_M2.png"), fig)
