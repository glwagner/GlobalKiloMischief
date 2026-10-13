#!/bin/bash
#SBATCH --job-name=gkm-download
#SBATCH --output=logs/download_data-%j.log
#SBATCH --account=bhcr-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00

# CDS returns several days of all variables in one file of about 7 GB, which is split into hourly files in memory:
# more than a login node allows. Compute nodes reach CDS.
#   sbatch slurm/download_data.sbatch [days]
source slurm/env.sh
julia --project=. download_data.jl "$@"
