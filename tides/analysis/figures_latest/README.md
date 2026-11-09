# Latest figures (tides SAL study, single-layer 0.36° model)

"calibrated" config: MOM6 aw depth (min 5 m, no cell removal), wave drag JSL_RH × 1.0 semi-implicit, flat biharmonic
A = Δ⁴/τ with τ = 1 d, ThicknessWeightedEnergy Coriolis (face h), Cd 0.003, g 9.8, M2 only, Δt 5 min, 60 days, fit days
30–60, SH n=40 PREM where SH. "native" = same numerics on our ETOPO depth (regridded on an 11000 m grid; layer = deepest cell).
Env: analysis/matrix/calibrated_env.sh (native: unset EXTERNAL_DEPTH). Run names below omit the _harmonics.jld2 suffix.

- matrix_amplitude.png: EOT20 M2 amplitude (top) and |Z − Z_EOT20| for {calibrated, native} × {none, scalar 0.1, SH n=40 PREM, convolution PREM}; scores in titles. Runs: analysis/matrix/runs_latest.tsv
  (barotropic_tides_0.36deg_dt5min_{none,scalar,spherical_harmonic}_wavedragjsl1_ThicknessWeightedEnergy_{cal,nat}_tau1_hf; convolution, SAL every step: barotropic_tides_0.36deg_dt5min_convolution_prem_wavedragjsl1_ThicknessWeightedEnergy_{cal,nat}_tau1_hf_acc_step, jobs 4158/4162).
- matrix_bars.png: errors (EOT20 deep/global, TPXO9 deep, gauges > 200 m) and cost (wall s per model day on a T4, SAL overhead) for the same runs.
- m2_overview.png: calibrated row: EOT20 with TICON-3 gauges, model vs gauge amplitudes, and |Z − Z_EOT20| for none/scalar/SH/convolution (every step). Gauge scores here use m2_overview.jl's own mask (our ETOPO wet cells), so they differ slightly from the matrix.
- kernel_summary.png: Eq. 16 ÷ 3 vs MOM6-coded kernel and their effect on deep error (historical numbers from the χ = 1.5 / MOM6 Table runs; unchanged, re-rendered).
- calibration_surface.png: (τ, χ) sweep on the calibrated base (Δt 5 min, SH): EOT20 deep (tuned), TPXO9 deep and gauges > 200 m (held out); red = chosen χ 1.0, τ 1 d (τ chosen for stability margin and conservatism, not on the gauges). Runs: barotropic_tides_0.36deg_dt5min_spherical_harmonic_wavedragjsl{χ}_ThicknessWeightedEnergy_cal_tau{τ}_hf.
- coriolis_energy_test.png: unforced inviscid energy E(t), native depth 5 m min, EnergyConserving vs ThicknessWeightedEnergy (logs etEnergyConserving-3988, etThicknessWeightedEnergy-3989).
- dt_convergence.png: deep amplitude ratio and rms difference vs the adjoint session's explicit model (level0_m036_t30_x0.8_nosal) for Δt 10/5/2.5 min (old Coriolis on both; runs barotropic_tides_0.36deg_dt{10,5,2.5min}_none_wavedragjsl0.8_EnergyConserving_cal_tau30).
- antarctic_mean_eta.png: fitted mean η south of 50°S before (default Coriolis) and after (ThicknessWeightedEnergy): runs barotropic_tides_0.36deg_dt10min_none_wavedrag1.5_s0_ref and ..._ThicknessWeightedEnergy_s0_twe.
