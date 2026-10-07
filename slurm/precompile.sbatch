#!/bin/bash
#SBATCH --job-name=gkso-precompile
#SBATCH --output=logs/precompile-%j.log
#SBATCH --account=bhcr-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00

# Precompile serially on the compute partition before any GPU or MPI job uses the environment.
export JULIA_NUM_PRECOMPILE_TASKS=16
source slurm/env.sh
julia --project=. -e 'using Pkg; Pkg.precompile(); using CUDA; CUDA.versioninfo(); using NumericalEarth, Oceananigans, ClimaSeaIce, MPI, NCCL, CopernicusMarine, CDSAPI'
