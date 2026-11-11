#!/bin/bash
#SBATCH --job-name=gkm-halobench
#SBATCH --output=logs/halo_benchmark-%j.log
#SBATCH --account=bhcr-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
#SBATCH --gres=gpu:4
#SBATCH --time=02:30:00

# NCCL against GPU-aware Cray MPI halo exchanges for the 1/12ᵒ configuration (bench/halo_benchmark.jl), one launch
# each, in the bench/ environment, whose MPI is the system Cray MPICH (the main environment's MPICH_jll is not
# GPU-aware).
#   sbatch bench/halo_benchmark.sbatch
source slurm/env.sh
source slurm/stage_inputs.sh 5
export LD_LIBRARY_PATH=/opt/cray/pe/mpich/9.0.1/ofi/gnu/12.3/lib:/opt/cray/pe/mpich/9.0.1/gtl/lib:/opt/cray/libfabric/2.3.1/lib64:/opt/cray/pe/pmi/default/lib
export MPICH_GPU_SUPPORT_ENABLED=1
julia --project=bench -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
for communication in nccl mpi; do
    srun -n 4 julia --project=bench bench/halo_benchmark.jl $communication
done
grep HALOBENCH logs/halo_benchmark-$SLURM_JOB_ID.log
