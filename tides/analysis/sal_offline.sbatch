#!/bin/bash
#SBATCH --partition=gpu-dev
#SBATCH --cpus-per-task=8
#SBATCH --output=/shared/home/greg/tides/logs/%x-%j.log
cd /shared/home/greg/tides
start=$(date +%s)
/shared/home/greg/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia -t 8 --project=sal_env analysis/sal_offline.jl \
    GlobalKiloMischief/tides/output/internal_tides_CUDA_Float64_0.125deg_Nz10_dt3min_tau30d.jld2 ${1:-10}
echo "OFFLINE EXIT $? after $(( $(date +%s) - start )) s"
