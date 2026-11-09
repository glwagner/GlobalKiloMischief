src = read(joinpath(@__DIR__,"check_aqua.jl"), String)
src = replace(src, "GPU()" => "CPU()", "days_run = 16.0" => "days_run = 0.05")
src = split(src, "a = zeros(Nx, Ny)")[1]
src = replace(src, "using NumericalEarth, Oceananigans, Oceananigans.Units, SpeedyTransforms, ConservativeRegridding, CUDA, Printf" => "using NumericalEarth, Oceananigans, Oceananigans.Units, Printf")
ARGS_ = nothing
empty!(ARGS); append!(ARGS, ["2","none","10"])
include_string(Main, src)
for n in 1:3
    time_step!(ocean.model, Δt)
    u = Array(interior(ocean.model.velocities.u)); v=Array(interior(ocean.model.velocities.v)); e=Array(interior(η))
    println(n," nan u ",count(isnan,u)," v ",count(isnan,v)," η ",count(isnan,e)," first ",findfirst(isnan,u)," ", findfirst(isnan,e), " maxu ", maximum(abs,filter(isfinite,u)))
end
Gu = Array(interior(ocean.model.timestepper.Gⁿ.u)); println("G nan ",count(isnan,Gu)," ",findfirst(isnan,Gu))
