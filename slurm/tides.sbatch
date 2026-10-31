#!/bin/bash
#SBATCH --job-name=gkm-tides
#SBATCH --output=logs/tides-%j.log
#SBATCH --account=bhcr-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00

# One GPU by default; for 1/12ᵒ on four: sbatch --ntasks=4 --gres=gpu:4 --mem=400G slurm/tides.sbatch 12 ...
#   sbatch slurm/tides.sbatch <cells per degree> <stop_days> <max_Δt_minutes> <cfl> <none | scalar | harmonic>
source slurm/env.sh
source slurm/stage_inputs.sh $2
# Precompile once before launching ranks, which would otherwise race on the shared depot.
julia --project=. -e 'using Pkg; Pkg.precompile()'
# One node: MPI (initialized even for one rank) talks over shared memory; it finds no network provider otherwise.
export FI_PROVIDER=shm
if [ "${SLURM_NTASKS:-1}" -gt 1 ]; then
    MPIEXEC=$(julia --project=. -e 'using MPI; print(MPI.mpiexec().exec[1])')
    $MPIEXEC -n $SLURM_NTASKS julia --project=. tides.jl "$@"
else
    julia --project=. tides.jl "$@"
fi
