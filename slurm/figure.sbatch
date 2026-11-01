#!/bin/bash
#SBATCH --job-name=gkm-figure
#SBATCH --output=logs/figure-%j.log
#SBATCH --account=bhcr-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --gres=gpu:1
#SBATCH --time=01:30:00

# Plotting packages live in their own environment.
#   sbatch slurm/figure.sbatch <run directory> <output.png> [stride]
source slurm/env.sh
env=/work/nvme/bhcr/glwagner/figure_env
julia --project=$env -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --project=$env figures/visualize_global_simulation_twelfth.jl "$@"
