#!/bin/bash
#SBATCH --partition=cpu-c6a24xlarge
#SBATCH --cpus-per-task=8
#SBATCH --output=/shared/home/greg/tides/logs/%x-%j.log
cd /shared/home/greg/tides
/shared/home/greg/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia -t 8 --project=sal_env analysis/m2_gauges.jl
echo "EXIT $?"
