#!/bin/bash
#SBATCH --job-name=gkm-cfl-tests
#SBATCH --output=logs/cfl_tests-%j.log
#SBATCH --account=bhcr-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
#SBATCH --gres=gpu:4
#SBATCH --time=02:00:00

# Four independent 1/6ᵒ runs side by side, one per GPU, each limited only by its CFL (the max_Δt cap is set high).
# Each run picks up from the latest checkpoint in its own directory, sixth_degree_dt<max_Δt>min_cfl<cfl>.
#   sbatch slurm/cfl_tests.sbatch <stop_days> <max_Δt_minutes> <cfl>...
source slurm/env.sh
source slurm/stage_inputs.sh $1
stop_days=$1; max_dt=$2; shift 2
for cfl in "$@"; do
    srun --ntasks=1 --gres=gpu:1 --cpus-per-task=16 --mem=95G --exact \
        --output=logs/cfl_tests-$SLURM_JOB_ID-cfl$cfl.log \
        julia --project=. sixth_degree.jl $stop_days $max_dt $cfl &
done
wait
