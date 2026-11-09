#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --time=4:00:00
#SBATCH --output=/shared/home/greg/tides/logs/%x-%j.log
cd /shared/home/greg/tides
/shared/home/greg/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia -t 4 --project=sal_env analysis/diag/sal_decomp.jl CUDA "$@"
echo "EXIT $?"
