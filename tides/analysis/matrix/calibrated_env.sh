# "calibrated" configuration (source before submitting): χ 1.0 JSL_RH, flat biharmonic τ 1 d, face-h thickness-weighted
# Coriolis, min 5 m, no cell removal, Δt 5 min (pass 5 to run.sbatch), M2 only, g 9.8, Cd 0.003, SH n=40 PREM when SH.
source /shared/home/greg/tides/analysis/matrix/calenv.sh
export WAVE_DRAG_COEFFICIENT=1.0 BIHARMONIC_TIMESCALE=1 SCALAR_LOADING_COEFFICIENT=0.1 CONVOLUTION_PRECISION=Float32 CONVOLUTION_ACCURATE=true
