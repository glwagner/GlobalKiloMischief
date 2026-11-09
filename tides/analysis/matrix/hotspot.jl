# Depth and M2 amplitude around the 20 m cell at (640, 50) (300.1ᵒE 62.3ᵒS) in a finished run
using JLD2, Printf
depth, λ, φ = load(joinpath(@__DIR__, "..", "diag", "grid_1000.jld2"), "depth", "λ", "φ")
out = joinpath(@__DIR__, "..", "..", "GlobalKiloMischief", "tides", "output")
for name in ARGS
    h = load(joinpath(out, name * "_harmonics.jld2"))
    Z = h["Z"]["M2"]
    println("== ", name, "  keys: ", keys(h))
    I, J = 636:644, 46:54
    println("depth (m), rows j = ", J)
    for j in reverse(J); println(j, ": ", join([@sprintf("%6.0f", depth[i, j]) for i in I], " ")); end
    println("|η| M2 (m)")
    for j in reverse(J); println(j, ": ", join([@sprintf("%6.2f", abs(Z[i, j])) for i in I], " ")); end
    if haskey(h, "Zu")
        Zu, Zv = h["Zu"]["M2"], h["Zv"]["M2"]
        println("|u| M2 (m/s) at u faces (i, j) = west face of cell i")
        for j in reverse(J); println(j, ": ", join([@sprintf("%6.2f", abs(Zu[i, j])) for i in I], " ")); end
        println("|v| M2 (m/s)")
        for j in reverse(J); println(j, ": ", join([@sprintf("%6.2f", abs(Zv[i, j])) for i in I], " ")); end
    end
    println("λ range ", extrema(λ[I, J]), " φ range ", extrema(φ[I, J]))
    # global ranking of the largest M2 amplitudes and depths
    o = sortperm(vec(abs.(Z)), rev = true)[1:10]
    for q in o; i, j = Tuple(CartesianIndices(Z)[q]); @printf("  |Z| %.2f m at (%d,%d) %.1fE %.1fN H %.0f\n", abs(Z[q]), i, j, λ[q], φ[q], depth[q]); end
end
