using JLD2, Printf
depth = load(joinpath(@__DIR__, "..", "diag", "grid_1000.jld2"), "depth")
h = load(joinpath(@__DIR__, "..", "..", "GlobalKiloMischief", "tides", "output", ARGS[1] * "_harmonics.jld2"))
for (I, J) in ((380:388, 10:16), (628:644, 46:54))
    for (lab, F) in (("depth", depth), ("mean η cm", 100h["mean"]), ("mean v cm/s", 100h["mean_v"]), ("mean u cm/s", 100h["mean_u"]))
        println(lab, " i=", I); for j in reverse(J); println(j, ": ", join([@sprintf("%6.0f", F[i, j]) for i in I], " ")); end
    end
end
