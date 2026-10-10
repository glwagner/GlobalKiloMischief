#!/bin/bash
#SBATCH --job-name=gkm-twelfth_degree
#SBATCH --output=logs/twelfth_degree-%j.log
#SBATCH --account=bhcr-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
#SBATCH --gres=gpu:4
#SBATCH --time=2-00:00:00

#   sbatch slurm/twelfth_degree.sbatch [stop_days] [max_Δt_minutes]
source slurm/env.sh
# Precompile once before launching ranks, which would otherwise race on the shared depot.
julia --project=. -e 'using Pkg; Pkg.precompile()'
# One node: MPI launches ranks and reduces scalars over shared memory, NCCL moves the GPU halos.
export FI_PROVIDER=shm
MPIEXEC=$(julia --project=. -e 'using MPI; print(MPI.mpiexec().exec[1])')
$MPIEXEC -n 4 julia --project=. twelfth_degree.jl "$@"
