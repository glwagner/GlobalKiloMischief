#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --time=3:00:00
#SBATCH --output=/shared/home/greg/tides/logs/%x-%j.log
cd /shared/home/greg/tides/GlobalKiloMischief/tides
echo "HOST $(hostname) SAL=$SELF_ATTRACTION_AND_LOADING CHI=$WAVE_DRAG_COEFFICIENT EPS=$MINIMUM_FRACTIONAL_CELL_HEIGHT TAG=$RUN_TAG OUT0=$OUTPUT_START_DAY CD=$BOTTOM_DRAG_COEFFICIENT"
start=$(date +%s)
/shared/home/greg/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia -t 4 --project=../../sal_env barotropic_tides.jl CUDA "$@"
echo "EXIT $? WALLTIME_SECONDS $(( $(date +%s) - start ))"
