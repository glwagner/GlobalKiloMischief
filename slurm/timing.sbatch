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
for variant in "$@"; do
    env=/work/nvme/bhcr/glwagner/sal_test/env
    [[ $variant == *+evaluated ]] && env=/work/nvme/bhcr/glwagner/sal_test/env_evaluated
    srun --ntasks=1 --gres=gpu:1 --cpus-per-task=16 --mem=95G --exact \
        --output=logs/timing-$SLURM_JOB_ID-$variant.log \
        julia --project=$env tides/timing.jl ${variant%+evaluated} &
    sleep 20  # launching all steps at once once failed with "Error configuring interconnect"
done
wait
