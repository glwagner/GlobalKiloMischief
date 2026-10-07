# Shared environment for GlobalKilometerScaleOceananigans jobs on DeltaAI (source this file).
export PATH=/u/glwagner/opt/julia-1.12.7/bin:$PATH
# Let CUDA.jl use its own runtime artifacts rather than the HPC SDK libraries on LD_LIBRARY_PATH.
unset LD_LIBRARY_PATH
export JULIA_NUM_THREADS=${JULIA_NUM_THREADS:-8}
export JULIA_NUM_PRECOMPILE_TASKS=${JULIA_NUM_PRECOMPILE_TASKS:-8}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-8}
export JULIA_PKG_PRECOMPILE_AUTO=0
# Downloaded and cached datasets (ERA5, GLORYS, JRA55, ETOPO) live on /work, not in the home quota.
export NUMERICALEARTH_DATA_DIRECTORY=/work/hdd/bhcr/glwagner/numericalearth_data
export RUN_DIRECTORY=${RUN_DIRECTORY:-/work/hdd/bhcr/glwagner/GlobalKilometerScaleOceananigans}
# Copernicus Marine credentials (COPERNICUSMARINE_SERVICE_USERNAME / _PASSWORD) for GLORYS.
[ -f "$HOME/.copernicusmarine.env" ] && source "$HOME/.copernicusmarine.env"
cd /u/glwagner/GlobalKilometerScaleOceananigans
mkdir -p logs "$RUN_DIRECTORY" "$NUMERICALEARTH_DATA_DIRECTORY"
echo "host: $(hostname)  job: ${SLURM_JOB_ID:-none}  start: $(date)"
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader 2>/dev/null || true
