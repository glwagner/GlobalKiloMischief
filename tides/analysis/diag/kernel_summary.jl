# Our Eq. 16 ÷ 3 kernel vs the kernel as coded in MOM6, and what each does to deep M2 error in both models.
# Scores are hand-entered from m2_atlases.jl and mom6_repro/score.txt.
using CairoMakie
n = [2, 3, 4, 6, 20, 40, 200.]; r = [.56, .77, .87, .92, .78, .68, .54]
eq = @. (7.21196 - 8.8 / n) / 3; prem = eq ./ r; code = @. (5.21196 + 8.865 / n) / 3
fig = Figure(size = (1250, 500), fontsize = 16)
a = Axis(fig[1, 1]; xscale = log10, xlabel = "degree n", ylabel = "kernel / PREM", title = "Two different \"paper\" kernels")
hlines!(a, 1; color = :gray)
scatterlines!(a, n, r; color = :dodgerblue4, label = "Our \"paper\" kernel: Eq. 16 ÷ 3")
scatterlines!(a, n, code ./ prem; color = :orangered3, marker = :rect, label = "As coded in MOM6: (5.21 + 8.87/n) / 3")
axislegend(a; framevisible = false)
labels = ["PREM-like\nSH n=40", "Eq. 16 ÷ 3\nSH n=40", "as coded\nSH n=40", "as coded\nconvolution"]
xo = [1, 2, 3]; o = [8.56, 7.65, 11.53]; xm = [1, 3, 4]; m = [8.11, 5.97, 5.96]
b = Axis(fig[1, 2]; xticks = (1:4, labels), ylabel = "deep M2 error (cm), area mean", title = "The coded kernel helps MOM6 but hurts ours")
barplot!(b, xo .- .2, o; width = .38, color = :dodgerblue4, label = "Oceananigans 0.36°, χ=1.5, vs EOT20", bar_labels = :y)
barplot!(b, xm .+ .2, m; width = .38, color = :goldenrod, label = "MOM6 p04arc9, Table χ, vs TPXO9", bar_labels = :y)
ylims!(b, 0, 14.5); axislegend(b; position = :rt, framevisible = false)
save(joinpath(@__DIR__, "kernel_summary.png"), fig)
