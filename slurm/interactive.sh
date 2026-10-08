#!/bin/bash
# Persistent Julia session on a GPU node for development, inside tmux so it survives disconnects:
#   tmux new -d -s gkm-repl 'bash slurm/interactive.sh [gpus=1] [hours=2]'
#   tmux attach -t gkm-repl
gpus=${1:-1}
hours=${2:-2}
srun --account=bhcr-dtai-gh --partition=ghx4-interactive --nodes=1 --ntasks=1 --gpus-per-node=$gpus \
     --cpus-per-task=16 --mem=200g --time=$(printf "%02d:00:00" $hours) --pty \
     bash -c 'source slurm/env.sh && julia --project=.'
