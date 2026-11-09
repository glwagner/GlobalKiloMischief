#!/bin/bash
#SBATCH --partition=cpu-c6a24xlarge
#SBATCH --cpus-per-task=96
#SBATCH --output=/shared/home/greg/tides/logs/%x-%j.log
# Figs. 1–3 of Chen et al. (2026) at 0.36ᵒ:  sbatch --export=ALL sal_convolution.sbatch <1/8ᵒ run> <day>
cd /shared/home/greg/tides
start=$(date +%s)
/shared/home/greg/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia -t 96 --project=sal_env analysis/sal_convolution.jl "$@"
echo "CONVOLUTION EXIT $? after $(( $(date +%s) - start )) s"
