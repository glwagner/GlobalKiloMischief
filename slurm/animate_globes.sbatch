#!/bin/bash
#SBATCH --job-name=gkm-globes
#SBATCH --output=logs/animate_globes-%j.log
#SBATCH --account=bhcr-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00

# The two-globe paper figure from the last snapshot, then an MP4 of every snapshot: CairoMakie draws on one core, so
# PROCESSES renderers (default 32) each take every PROCESSES-th frame, and ffmpeg joins the frames at FRAMERATE (24).
#   sbatch slurm/animate_globes.sbatch <run directory> <output prefix>   # writes <prefix>.png and <prefix>.mp4
source slurm/env.sh
env=/work/nvme/bhcr/glwagner/figure_env
run=$1
prefix=$2
processes=${PROCESSES:-32}
framerate=${FRAMERATE:-24}
frames=$prefix-frames
julia --project=$env -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

julia --project=$env figures/visualize_global_simulation_twelfth.jl $run $prefix.png 2

count=$(julia --project=$env figures/animate_global_simulation_twelfth.jl $run count | tail -1)
echo "rendering $count frames with $processes processes"
for p in $(seq 1 $processes); do
    JULIA_NUM_THREADS=1 julia --project=$env figures/animate_global_simulation_twelfth.jl $run $frames $p $processes 2 \
        > logs/animate_globes-$SLURM_JOB_ID-$p.log 2>&1 &
done
wait
echo "rendered $(ls $frames | grep -c '^frame_') of $count frames"

julia --project=$env -e "using FFMPEG_jll; run(\`\$(ffmpeg()) -y -loglevel error -framerate $framerate -i $frames/frame_%04d.png -c:v libx264 -pix_fmt yuv420p -crf 18 $prefix.mp4\`)"
ls -la $prefix.png $prefix.mp4
