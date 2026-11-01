#!/bin/bash
#SBATCH --job-name=gkm-tides-timing
#SBATCH --output=logs/timing-%j.log
#SBATCH --account=bhcr-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
#SBATCH --gres=gpu:4
#SBATCH --time=01:30:00

# Up to four tides/timing.jl configurations side by side, one per GPU. The environment develops NumericalEarth's
# glw/distributed-harmonic-sal, whose SAL update is split into the functions the timing calls.
#   sbatch slurm/timing.sbatch <configuration>...
source slurm/env.sh
source slurm/stage_inputs.sh 3
export FI_PROVIDER=shm
env=/work/nvme/bhcr/glwagner/sal_test/env
julia --project=$env -e 'using Pkg; Pkg.precompile()'
for variant in "$@"; do
    srun --ntasks=1 --gres=gpu:1 --cpus-per-task=16 --mem=95G --exact \
        --output=logs/timing-$SLURM_JOB_ID-$variant.log \
        julia --project=$env tides/timing.jl $variant &
    sleep 20  # launching all steps at once once failed with "Error configuring interconnect"
done
wait
