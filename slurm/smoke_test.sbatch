#!/bin/bash
#SBATCH --job-name=gos-smoke_test
#SBATCH --output=logs/smoke_test-%j.log
#SBATCH --account=bhcr-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00

#   sbatch slurm/smoke_test.sbatch [args...]
source slurm/env.sh
julia --project=. smoke_test.jl "$@"
