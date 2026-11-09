using NumericalEarth, Oceananigans, JLD2, Printf
Nx, Ny = 1000, 472
grid = TripolarGrid(CPU(); size = (Nx, Ny, 1), halo = (7, 7, 4), z = (-6000, 0))
H = -interior(regrid_bathymetry(grid; minimum_depth = 20, interpolation_passes = 3, major_basins = 1), :, :, 1)
h = load(joinpath(@__DIR__, "..", "..", "GlobalKiloMischief", "tides", "output", ARGS[1] * "_harmonics.jld2"))
η = h["η_snapshots"]["day010"]
for (ic, jc) in ((260, 167), (633, 279))
    println("η (m) around ($ic, $jc), rows j descending"); for j in jc+3:-1:jc-3; println(join([@sprintf("%7.2f", η[i, j]) for i in ic-3:ic+3], " ")); end
    println("H (m)"); for j in jc+3:-1:jc-3; println(join([@sprintf("%7.0f", H[i, j]) for i in ic-3:ic+3], " ")); end
end
