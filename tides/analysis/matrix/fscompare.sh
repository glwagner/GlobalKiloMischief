#!/bin/bash
O=/shared/home/greg/tides/GlobalKiloMischief/tides/output
cd /shared/home/greg/tides/egbert_ray/solve
for n in dt10min_none_wavedragjsl0.8_ThicknessWeightedEnergy_cal_tau30 dt5min_none_wavedragjsl0.8_ThicknessWeightedEnergy_cal_tau30 dt2.5min_none_wavedragjsl0.8_ThicknessWeightedEnergy_cal_tau30 dt10min_none_wavedragjsl0.8_ThicknessWeightedEnergy_cal_tau30_fscosine dt10min_none_wavedragjsl0.8_ThicknessWeightedEnergy_cal_tau30_fsconstant dt10min_none_wavedragjsl0.8_ThicknessWeightedEnergy_cal_tau30_cfl0.35; do echo "=== $n"; /shared/home/greg/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia -t 8 --project=/shared/home/greg/tides/sal_env tp036_compare.jl m036_t30_x0.8_nosal $O/barotropic_tides_0.36deg_${n}_harmonics.jld2 2>&1 | grep -E "^deep H|ratio"; done
