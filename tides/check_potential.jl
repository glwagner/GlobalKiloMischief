# Compare EquilibriumTidalPotential (M2, S2, K1, O1) with Tidejinks' full astronomical potential at a few points over
# 30 days of 2015: correlation and amplitude ratio (a sign or phase error shows up as a low or negative correlation).
#   julia --project=. tides/check_potential.jl

using Oceananigans
using Oceananigans.Units
using Dates: Dates, DateTime, Day
using Printf
using Statistics
using SPICE: SPICE
using Tidejinks: Tidejinks
using NumericalEarth

const start_date = DateTime(2015, 1, 1)
const Adapt = NumericalEarth.Oceans.Adapt
include(joinpath(@__DIR__, "stratified_tides.jl"))

grid = LatitudeLongitudeGrid(CPU(); size = (4, 4, 1), longitude = (0, 360), latitude = (-60, 60), z = (-1, 0))
dates = (start_date, start_date + Day(31))
cd(mkpath(joinpath(ENV["NUMERICALEARTH_DATA_DIRECTORY"], "SPICE"))) do
    isfile("kernels.txt") || Tidejinks.wrangle_spice_kernels("kernels.txt", "kernels")
    SPICE.furnsh("kernels.txt")
end
reference = Tidejinks.TidalPotential(grid, dates; reference_frame = "IAU_EARTH")
Φ = EquilibriumTidalPotential(start_date; ramp = 1)  # ramped up long before day 3

times = 3days:1hour:30days
for (λ, φ) in ((0, 0), (90, 30), (200, -45), (300, 60))
    a = [Φ(λ, φ, t) for t in times]
    b = [reference(λ, φ, t) for t in times]
    @printf("λ = %3d, φ = %3d: correlation %.3f, std ratio (4 constituents / Tidejinks) %.3f\n", λ, φ, cor(a, b), std(a) / std(b))
end

