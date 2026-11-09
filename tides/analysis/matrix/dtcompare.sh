#!/bin/bash
O=/shared/home/greg/tides/GlobalKiloMischief/tides/output
cd /shared/home/greg/tides/egbert_ray/solve
for dt in 10 5 2.5; do echo "=== dt $dt"; /shared/home/greg/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia -t 8 --project=/shared/home/greg/tides/sal_env tp036_compare.jl m036_t30_x0.8_nosal $O/barotropic_tides_0.36deg_dt${dt}min_none_wavedragjsl0.8_EnergyConserving_cal_tau30_harmonics.jld2 2>&1 | grep -v Warn | grep -E "deep|global|rms|ratio"; done
