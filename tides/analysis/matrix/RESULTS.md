# SAL matrix: {MOM6-matched, native} × {none, scalar 0.1, SH n=40 PREM, convolution PREM}

Single-layer global tide model `GlobalKiloMischief/tides/barotropic_tides.jl`, 0.36ᵒ TripolarGrid (1000 × 472), 60 days,
online harmonic fit over days 30–60. Errors are rms|ΔZ|/√2 in cm, plain area mean (`matrix/score_matrix.jl`), on cells
wet in both our bathymetry and the regridded MOM6 one. Deep = depth > 1000 m and |φ| < 66ᵒ. "amp" = rms(|m|−|o|)/√2,
"phase" = sqrt(total² − amp²), "offset" = area-weighted mean phase difference model − atlas (ᵒ). Gauges: TICON-3
coastal M2, nearest wet cell within 1.5 cells (2137 gauges; 349 at model depth > 200 m). Cost on one NVIDIA T4
(gpu-dev, g4dn.2xlarge).

Figures: `matrix_amplitude.png` (M2 amplitude, errors in titles), `matrix_bars.png` (errors, cost).

## Matrix

| config | SAL | EOT20 deep | amp / phase | offset ᵒ | EOT20 global | TPXO9 deep | TPXO9 global | gauges all | gauges >200 m | wall s / model day | ms / step | SAL overhead s / day |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| MOM6-matched | none | 24.36 | 8.33 / 22.89 | +45.2 | 26.14 | 24.21 | 26.30 | 66.1 | 34.5 | 31.2 | 65 | – |
| MOM6-matched | scalar 0.1 | 8.50 | 5.27 / 6.67 | −5.7 | 14.04 | 8.73 | 15.49 | 56.3 | 22.4 | 31.4 | 65 | 0.2 |
| MOM6-matched | SH n=40 PREM | **6.04** | 3.76 / 4.73 | −3.5 | 9.86 | 6.30 | 11.54 | 49.4 | 18.6 | 54.0 | 113 | 22.8 |
| MOM6-matched | convolution PREM ¹ | 6.11 | 4.07 / 4.56 | −1.2 | 9.85 | 6.29 | 11.32 | 48.7 | 18.2 | 190.8 | 398 | 159.6 |
| native | none | 17.52 | 5.55 / 16.62 | +34.7 | 19.24 | 17.39 | 19.65 | 53.5 | 26.7 | 12.0 | 83 | – |
| native | scalar 0.1 | 10.52 | 6.25 / 8.47 | −8.6 | 16.43 | 10.69 | 17.57 | 55.4 | 21.6 | 12.1 | 84 | 0.1 |
| native | SH n=40 PREM | 8.31 | 5.96 / 5.78 | −5.3 | 12.65 | 8.48 | 13.98 | 44.7 | 19.3 | 18.6 | 129 | 6.6 |
| native | convolution PREM ¹ | 9.01 | 6.82 / 5.88 | −4.1 | 12.94 | 9.13 | 14.08 | 44.6 | 19.4 | 165.2 | 1147 | 153.2 |

¹ Float32 with `CONVOLUTION_ACCURATE=true` (fold duplicates removed, uniform-response subtraction; transfer function
identical to Float64 to 4 digits, see below), SAL recomputed every 30 min and held in between (`SAL_UPDATE_INTERVAL`
= 10 steps of 3 min / 3 steps of 10 min). Fair comparison with SH held the same way:

| config | SH n=40 held 30 min | EOT20 deep | offset ᵒ | TPXO9 deep | gauges all / >200 |
|---|---|---|---|---|---|
| MOM6-matched | SH, held 10 steps | 6.02 | −0.7 | 6.19 | 48.4 / 18.1 |
| native | SH, held 3 steps | 8.84 | −3.6 | 8.95 | 43.8 / 19.4 |

Holding the SAL for 30 min shifts the mean phase by ~2–3ᵒ in both rows; the deep error is unchanged in the MOM6-matched
row (6.04 → 6.02) but worsens by 0.5 cm in the native row (8.31 → 8.84). Like for like (both held), convolution vs SH
is 6.11 vs 6.02 (MOM6-matched) and 9.01 vs 8.84 (native): **the exact PREM convolution is no better than SH n=40 at this
resolution, at ~25–50× the SAL cost**, and its smaller phase offset in the main table is mostly the hold lag.

Run names: `runs.tsv`; full scores: `scores.tsv`; costs: `costs.tsv`.

### SAL cost (synchronised timing of `compute_self_attraction_and_loading!` alone, 20 calls)

| method | ms per update | notes |
|---|---|---|
| scalar | 0 | folded into g; whole-run difference vs none ≤ 0.2 s/day (noise) |
| SH n=40 (48 × 96 Gauss grid, conservative regridding) | 21–22 | called at every time-stepper substage (~2 per step), so +46 ms/step |
| convolution Float32 (default) | 3,070 | direct O(N²), N ≈ 3.0 × 10⁵ wet cells |
| convolution Float32 accurate | 3,230–3,360 | +8% for the subtraction |
| convolution Float64 | 50,000–57,000 in the model; 130 s on the all-wet aquaplanet | T4 FP64 is 1/32 rate: infeasible (60 days with an update every step at Δt 10 min ≈ 140 days of wall time) |

## Configurations

MOM6-matched (`matrix/m6env.sh`), from Chen et al.'s p04arc9 MOM_input / MOM_override / MOM_parameter_doc.all:
- depth: `depth_p04arc9_aw_moreocean_ocnpct0p5.nc` regridded by Voronoi binning (`diag/regrid_mom6.jl`: each MOM6 h-point
  to our nearest cell centre, wet if ≥ half its points are wet (MOM6 MASKING_DEPTH = −5 m), mean over wet points; 1.59
  points per cell; cells with none take the nearest MOM6 point within 0.5ᵒ). 301,604 wet cells (ours 294,458).
  Layer as deep as the deepest cell (9314 m; ours clips at 6000 m). No isolated-cell removal.
- minimum depth 5 m: stable (1005 cells deepened). MOM6 itself keeps 744 cells with bottom up to 5 m above sea level
  wet (thickness D + η); we cannot, so those are 5 m deep.
- wave drag: JSL_RH × BT_WAVE_DRAG_SCALE 0.8 (the SH run's value; their convolution runs used 0.7), r = 0.8 JSL_RH / H,
  semi-implicit. Max r = 1.02 × 10⁻⁵ s⁻¹ (1/r = 27 h; the file is pre-limited to 1/24 h at scale 1): r Δt = 0.0018 at
  Δt = 180 s, so explicit would have been stable too; implicit changes nothing (also checked in our config: 8.55 vs 8.56).
- viscosity: Laplacian KH = 0.03 Δ, biharmonic AH = max(0.01 Δ³, 0.06 Δ⁴ |S|), Δ² = 2Δx²Δy²/(Δx²+Δy²), bounded
  Δ²/(8Δt), Δ⁴/(64Δt) (RESOLN_SCALED_KH is inert with one layer). KH ≤ 1234 m² s⁻¹, AH ≤ 7 × 10¹¹ m⁴ s⁻¹.
- bottom drag Cd 0.003 (with one layer HBBL only sets the velocity sampled, so same as ours), g = 9.8, Δt = 180 s,
  SADOURNY75_ENERGY → EnergyConserving Coriolis, M2-only equilibrium forcing (TIDE_M2), as in MOM6.
- not matched: MOM6's MAXVEL = 6 m/s clipping, BOUND_CORIOLIS, BT_NONLIN_STRESS details, split-RK2 time stepping.

Native, our own Oceananigans setup (`matrix/bestenv.sh`): our ETOPO bathymetry (min 20 m, isolated cells removed), wave drag JSL_RH × 1.1 on our
depth, uniform Laplacian ν = 1000 m² s⁻¹, Cd 0.003, g = 9.80665, Δt = 10 min, full Tidejinks potential (all constituents).

### Sweep (SH n=40 PREM, our bathymetry, 14 runs; tuned on EOT20 deep only)

| wave drag | viscosity | EOT20 deep | TPXO9 deep (held out) | gauges all / >200 (held out) |
|---|---|---|---|---|
| Egbert–Ray C χ = 1 | ν 1000 | 9.33 | 9.57 | 47.7 / 20.1 |
| χ = 1.5 | ν 1000 | 8.56 | 8.76 | 45.4 / 19.1 |
| χ = 1.5, implicit | ν 1000 | 8.55 | 8.75 | 45.5 / 19.1 |
| χ = 2 | ν 1000 | 8.63 | 8.78 | 44.0 / 18.7 |
| χ = 3 | ν 1000 | 9.89 | 9.96 | 43.0 / 18.7 |
| χ = 1.5 | ν 250 | 8.55 | 8.75 | 45.4 / 19.1 |
| χ = 1.5 | ν 500 | 8.55 | 8.75 | 45.4 / 19.1 |
| χ = 1.5 | MOM6 grid-scaled | 8.49 | 8.69 | 45.2 / 19.1 |
| χ = 2 | MOM6 grid-scaled | 8.56 | 8.72 | 43.9 / 18.6 |
| JSL_RH × 0.6 | ν 1000 | 8.95 | 9.19 | 47.9 / 20.3 |
| JSL_RH × 0.8 | ν 1000 | 8.35 | 8.56 | 46.3 / 19.7 |
| JSL_RH × 0.8 | MOM6 grid-scaled | 8.31 | 8.53 | 46.1 / 19.7 |
| **JSL_RH × 1.1** | **ν 1000** | **8.31** | **8.48** | **44.7 / 19.3** |
| JSL_RH × 1.6 | ν 1000 | 9.43 | 9.52 | 43.8 / 19.4 |

Justification. Viscosity is irrelevant at this resolution for M2 (250–1000 m² s⁻¹ identical to 0.01 cm): the tide is
not viscously controlled, so we keep 1000 m² s⁻¹, the smallest value that is cheap insurance against grid noise
next to steep partial cells. Wave drag controls the error: the deep error has a broad minimum where the conversion
to internal tides is ~ the Jayne–St. Laurent estimate itself. JSL_RH (computed on 1/12ᵒ topography) at scale ≈ 1
(1.1) is the unscaled physical value and beats our own 0.5ᵒ Egbert–Ray C field at its tuned χ = 1.5 (which compensates
for unresolved roughness on coarse tiles); the held-out TPXO9 and gauge scores rank the same way (TPXO9 8.48 vs
8.76; gauges 44.7 vs 45.4). Larger drag improves gauges (shelf amplitudes) but degrades the deep ocean (amplitude
part grows), the usual trade-off.

**Caveat: the MOM6-matched row beats our native setup (6.04 vs 8.31 deep EOT20).** The tuned sweep was all on our
bathymetry; the differences between rows are bathymetry (MOM6 aw depth, 5 m minimum, deeper trenches, no isolated
cell removal), Δt 3 vs 10 min, M2-only forcing, the MOM6 viscosity, g, and the Coriolis scheme. Not yet separated
(see "What remains"); on gauges the native row is better (44.7 vs 49.4 all, similar > 200 m).

## Coded paper kernel in the MOM6-matched configuration

SH n=40 with `LOADING_LOVE_NUMBERS=paper_code` (Chen et al.'s coded kernel incl. ÷3): EOT20 deep **8.65** (amp 5.00,
phase 7.06, offset −12.6ᵒ) vs 6.04 with PREM; TPXO9 8.86 vs 6.30. So the coded kernel hurts in our model even with
MOM6's inputs (as in our own config: 11.5 vs 8.56), unlike in MOM6 (5.97 vs 8.11). Matching depth, drag, viscosity,
min depth and g does not explain MOM6's preference for it; the remaining suspects are in MOM6's SAL implementation
or dynamics, not the inputs.

## Defects found (with evidence)

1. **Oceananigans' biharmonic closure + true-depth partial cells is unstable** (MOM6 viscosity via
   HorizontalScalarBiharmonicDiffusivity): η of 6000 m in deep water after 0.02 days with 5 m cells (job 3793), NaN at
   day ~2 with our 20 m bathymetry (3804). Both Laplacians are thickness-weighted, so a thin cell next to a deep one
   amplifies the inner Laplacian by the depth ratio. Fixed by applying the MOM6-style viscosity as a forcing with flat
   (not thickness-weighted) 2D operators, free slip, `VISCOSITY=mom6|flat` (the header's note that the biharmonic
   "drives steady grid-scale η of tens to hundreds of metres" is probably the same thing).
2. **Convolution: tripolar fold row counted twice.** Cells (i, 472) and (1001−i, 472) coincide (d ≈ 5 × 10⁻¹⁷ rad;
   310 wet pairs in our mask, 154–155 in the run grids; `check_duplicates.jl`). They are double sources and, in
   Float64, interact at 1 − cos γ = eps → dG ~ 10²³ (the "rel residual 4.5" of the old Float64 aquaplanet test).
   `remove_duplicates=true` zeroes the duplicate's source weight and skips coincident pairs.
3. **Convolution Float32 low-degree error is cancellation, not accumulation.** Transfer at n = 1, 2, 3: Float32
   0.895, 0.954, 0.989; Float64 0.979, 1.000, 1.000. `subtract_uniform=true` (sum Ωₘ(pₘ − pₙ)(xₘ − xₙ) dG in Float32,
   plus pₙ × the uniform-load response precomputed once in Float64) gives 0.979, 1.000, 1.000 — identical to Float64
   at all tested degrees (1, 2, 3, 5, 10, 40) — for +8% cost instead of ×17 (`conv_transfer.jl`, job 3820). Adjoint
   check of the unchanged transpose: 3–5 × 10⁻³ (Float32 roundoff of random fields), 2.7 × 10⁻³ in Float64.
   The residual ~0.087 at all degrees is the same for all variants (metric floor of the comparison).
   Also: with `components = :grid` four aquaplanet cells at the tripolar poles have NaN e¹ (land in real masks).
4. **SAL_UPDATE_INTERVAL via an UpdateStateCallsite + IterationInterval does not hold the SAL** (cost and result
   unchanged); it now uses an end-of-step callback. Also, the SH SAL callback runs at every time-stepper substage
   (~2 per step), doubling its cost.
5. **Large tidally rectified mean flows near Antarctica** (open): the fitted mean η reaches 6.5 m (South Shetlands,
   300ᵒE 62ᵒS, and Ross/Amundsen coast 205–215ᵒE 75ᵒS) with steady along-slope currents up to 12 m/s in 3000 m
   water; > 0.5 m/s in ~2700 cells. This is the "20 m cell at 300ᵒE 62ᵒS reaching 6 m/s": its M2 amplitude is only
   0.39 m; the 6.5 m is the mean. Evidence (`meanflow.jl`, `peek.jl`, 30-day diagnostic runs): absent with zero
   forcing, 4.3 m at half M2 amplitude, 6.3 m at full (strongly sub-quadratic, saturating), the same with no SAL,
   M2-only forcing, and EnergyConserving / ActiveWeightedEnstrophy Coriolis; not the partial-cell hydrostatic
   pressure of the uniform buoyancy (no correlation of mean η with H). Not root-caused; it also appears in the
   MOM6-matched runs (max|η| ~ 9–12 m at 212ᵒE 75.7ᵒS). Mostly outside the deep |φ| < 66ᵒ mask, but it adds drag
   locally and the Weddell/Filchner M2 amplitude (> 1.5 m, bottom row of the figure) is suspect.
6. **Arctic resonance:** M2 amplitude 5–6 m in the Lincoln Sea / north of Greenland (81ᵒN, 290–295ᵒE, near the
   tripolar fold) in our bathymetry (real tide ~0.3 m).
7. Δt in run.sbatch's documented call is 10 min (not 1 min); MOM6 uses 180 s. Wall times above are per model day.

## Infrastructure added (env options; defaults reproduce earlier results)

`barotropic_tides.jl`: `EXTERNAL_DEPTH=mom6`, `MINIMUM_DEPTH` (20), `GRAVITATIONAL_ACCELERATION` (9.80665),
`WAVE_DRAG_FIELD=jsl_rh`, `WAVE_DRAG_IMPLICIT`, `VISCOSITY=mom6|flat` (+ `MOM6_VISCOSITY_SCALE`),
`SAL_UPDATE_INTERVAL`, `CONVOLUTION_ACCURATE`, `M2_AMPLITUDE_SCALE` (M2_ONLY only); COST / SALCOST log lines; fit
files also store mean u, v. `MINIMUM_FRACTIONAL_CELL_HEIGHT` default becomes min(1e-3, MINIMUM_DEPTH / 2 layer depth)
(= 1e-3 for the defaults).
NumericalEarth-sal `ConvolutionLoadingOperator(...; remove_duplicates = false, subtract_uniform = false)`: new keywords,
defaults unchanged; `love_numbers` and `loading_gradient_transpose!` API unchanged (the transpose honours the
duplicate weights).

## What remains

- Separate why MOM6-matched beats native (bathymetry vs Δt vs M2-only vs min depth): 3–4 single-factor SH runs.
- Root-cause the rectified Antarctic mean flows (defect 5) and the Arctic resonance (6).
- Why MOM6 prefers the coded kernel: compare the SAL fields of both models for the same η (sal_decomp.jl part B on a
  MOM6 snapshot).
- Convolution updated every step (needs ~8 h per run at Δt 10 min on a T4, ~20 h at 3 min, or an H100) to remove the
  hold lag from the convolution column.
