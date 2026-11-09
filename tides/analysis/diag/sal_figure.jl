using CairoMakie
L = Dict(n=>v for (n,v) in ((1,-0.412),(2,.556),(3,.769),(4,.871),(5,.917),(8,.923),(10,.897),(15,.829),(20,.776),(30,.710),(40,.675),(80,.617),(200,.536)))
ns = sort(collect(keys(L)))
pw = [(1,.0007),(2,.0555),(3,.0501),(4,.1381),(5,.1925),(6,.095),(7,.053),(8,.0524),(9,.0378),(10,.0219),(20,.0064),(40,.0032)]
f = Figure(size=(900,380))
a = Axis(f[1,1], xscale=log10, xlabel="degree n", ylabel="(fit/3) / PREM kernel", title="Paper kernel vs PREM, with M2 power fraction")
hlines!(a,[1.0],color=:gray); scatterlines!(a,ns,[L[n] for n in ns],label="ratio"); barplot!(a,first.(pw),last.(pw),color=(:orange,.6),width=[.25n for (n,_) in pw],label="M2 η power"); axislegend(a,position=:rb)
b = Axis(f[1,2], xlabel="χ", ylabel="deep M2 error vs EOT20 (cm)", title="χ sweep, 0.36°")
scatterlines!(b,[.5,1,1.5,2,3],[11.31,9.34,8.56,8.63,9.89],label="SH n=40 PREM-like")
scatterlines!(b,[1.5,2,3],[7.65,7.78,9.31],label="SH n=40 paper kernel")
scatter!(b,[2],[8.87],label="convolution PREM",marker=:utriangle,markersize=14); scatter!(b,[2],[7.80],label="convolution paper",marker=:diamond,markersize=14)
axislegend(b,position=:ct); save(joinpath(@__DIR__,"sal_attribution.png"),f)
