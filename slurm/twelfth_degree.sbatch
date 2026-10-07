#!/bin/bash
#SBATCH --job-name=gkso-twelfth_degree
#SBATCH --output=logs/twelfth_degree-%j.log
#SBATCH --account=bhcr-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
#SBATCH --gres=gpu:4
#SBATCH --time=2-00:00:00

#   sbatch slurm/twelfth_degree.sbatch [stop_days]
# The environment must already be precompiled (slurm/precompile.sbatch): ranks precompiling concurrently race.
source slurm/env.sh
MPIEXEC=$(julia --project=. -e 'using MPI; print(MPI.mpiexec().exec[1])')
$MPIEXEC -n 4 julia --project=. twelfth_degree.jl "$@"
