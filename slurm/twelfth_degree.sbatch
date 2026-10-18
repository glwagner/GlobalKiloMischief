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
source slurm/stage_inputs.sh ${1:-30}  # stop_days
# Precompile once before launching ranks, which would otherwise race on the shared depot.
julia --project=. -e 'using Pkg; Pkg.precompile()'
# One node: MPI launches ranks and reduces scalars over shared memory, NCCL moves the GPU halos.
export FI_PROVIDER=shm
MPIEXEC=$(julia --project=. -e 'using MPI; print(MPI.mpiexec().exec[1])')

# A 4-rank run hung after 10 steps without an error (job 3354360). If the log stops growing for 15 minutes once
# stepping has begun, ask every rank for a profile of where its threads are (Julia prints one on SIGUSR1).
log=logs/twelfth_degree-$SLURM_JOB_ID.log
( while sleep 300; do
      grep -q "initial time step complete" $log || continue
      [ $(( $(date +%s) - $(stat -c %Y $log) )) -gt 900 ] && pkill -USR1 -u $USER -f "twelfth_degree.jl" && sleep 900
  done ) &

$MPIEXEC -n 4 julia --project=. twelfth_degree.jl "$@"
