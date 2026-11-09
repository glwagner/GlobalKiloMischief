# Equilibrium M2 from Tidejinks TidalPotential at a few points, converted with barotropic_m2.jl's Greenwich convention.
using Oceananigans, Tidejinks, Dates, LinearAlgebra, Printf
import SPICE
cd(joinpath(@__DIR__, "..", "..", "GlobalKiloMischief", "tides"))
SPICE.furnsh("kernels.txt")
start_date = DateTime(2015, 1, 1)
grid = RectilinearGrid(size=(1,1,1), extent=(1,1,1))
Φ = Tidejinks.TidalPotential(grid, (start_date, start_date + Day(31)); reference_frame = "IAU_EARTH")
first_day, last_day = 5, 20
periods = (12.4206012, 12.0, 23.9344696, 25.8193417, 12.65834751, 11.96723606)
ω = [2π / (3600T) for T in periods]
t = collect(first_day*86400:3600:last_day*86400)
A = hcat(ones(length(t)), [cos(w*x) for x in t, w in ω], [sin(w*x) for x in t, w in ω]); P = pinv(A); K = length(ω)
julian_centuries(date) = (Dates.datetime2julian(date) - 2451545.0) / 36525
T₀ = julian_centuries(start_date)
V₀ = mod(2(280.4664567 + 36000.76983T₀) - 2(218.3164477 + 481267.88123421T₀), 360)
N = 125.04452 - 1934.136261julian_centuries(start_date + Hour(300))
u = -2.14sind(N); f = 1.0004 - 0.0373cosd(N) + 0.0002cosd(2N)
g = 9.80665
for (λ, φ) in ((0.0, 0.0), (90.0, 0.0), (0.0, 45.0), (-150.0, 20.0))
    η = [-Φ(λ, φ, x) / g for x in t]
    c = P * η
    Z = c[2] - im * c[2+K]
    G = mod(-rad2deg(angle(Z)) + V₀ + u, 360)
    Aeq = abs(Z) / f
    theory = 0.2423 * 0.69 * cosd(φ)^2  # Cartwright M2 equilibrium amplitude × (1+k-h)
    @printf("λ=%6.1f φ=%5.1f: A=%.4f m (theory %.4f) G=%7.2f (expect %7.2f)\n", λ, φ, Aeq, theory, G, mod(-2λ, 360))
end
