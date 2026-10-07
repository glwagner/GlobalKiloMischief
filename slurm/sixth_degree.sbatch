#!/bin/bash
#SBATCH --job-name=gkm-sixth_degree
#SBATCH --output=logs/sixth_degree-%j.log
#SBATCH --account=bhcr-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --gres=gpu:1
#SBATCH --time=2-00:00:00

#   sbatch slurm/sixth_degree.sbatch [args...]
source slurm/env.sh
julia --project=. sixth_degree.jl "$@"
