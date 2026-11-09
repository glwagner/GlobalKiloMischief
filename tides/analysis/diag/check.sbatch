#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --time=6:00:00
#SBATCH --output=/shared/home/greg/tides/logs/%x-%j.log
cd /shared/home/greg/tides
J=/shared/home/greg/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia

for c in "2 none 10" "2 none 5" "2 scalar 10" "2 sh 10" "2 sh 5" "2 conv 10" "2 conv 5" "4 none 10" "4 sh 10" "4 conv 10"; do
 $J -t 4 --project=sal_env analysis/diag/check_aqua.jl $c 2>&1 | grep -E 'RESULT|rror' | cut -c1-200
done
echo EXIT
