#!/bin/bash
#SBATCH --job-name=gkm-tides-timing
#SBATCH --output=logs/timing-%j.log
#SBATCH --account=bhcr-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
#SBATCH --gres=gpu:4
#SBATCH --time=01:30:00

# Up to four tides/timing.jl configurations side by side, one per GPU. The environment develops NumericalEarth's
# glw/distributed-harmonic-sal, whose SAL update is split into the functions the timing calls; a configuration named
# with the suffix +evaluated uses glw/evaluated-barotropic-potential instead, which evaluates the tidal potential into
# a 2D field once per state update rather than inside the momentum tendency kernels.
#   sbatch slurm/timing.sbatch <configuration[+evaluated]>...
source slurm/env.sh
source slurm/stage_inputs.sh 3
export FI_PROVIDER=shm
for env in /work/nvme/bhcr/glwagner/sal_test/env /work/nvme/bhcr/glwagner/sal_test/env_evaluated; do
    julia --project=$env -e 'using Pkg; Pkg.precompile()'
done
# Each configuration is a plain process on its own GPU: concurrent srun steps fail to launch on DeltaAI
# ("Error configuring interconnect"), even 20 s apart.
gpu=0
for variant in "$@"; do
    env=/work/nvme/bhcr/glwagner/sal_test/env
    [[ $variant == *+evaluated ]] && env=/work/nvme/bhcr/glwagner/sal_test/env_evaluated
    CUDA_VISIBLE_DEVICES=$gpu JULIA_NUM_THREADS=16 julia --project=$env tides/timing.jl ${variant%+evaluated} \
        > logs/timing-$SLURM_JOB_ID-$variant.log 2>&1 &
    gpu=$((gpu + 1))
done
wait
