#!/bin/bash
#SBATCH --job-name=gkm-animate
#SBATCH --output=logs/animate-%j.log
#SBATCH --account=bhcr-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00

# Plotting packages live in their own environment, so the simulation environment does not carry CairoMakie.
#   sbatch slurm/animate.sbatch <output.mp4> <run directory>...
source slurm/env.sh
env=/work/nvme/bhcr/glwagner/animation/env
julia --project=$env -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --project=$env animate_surface_speed.jl "$@"
