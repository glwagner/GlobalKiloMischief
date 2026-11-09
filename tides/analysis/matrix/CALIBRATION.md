# Calibration of the single-layer 0.36ᵒ tide model (step 0–3)

Model: `GlobalKiloMischief/tides/barotropic_tides.jl`. Scores: `matrix/score_matrix.jl` (rms|ΔZ|/√2, cm, plain area
mean; deep = H > 1000 m, |φ| < 66ᵒ; mask wet in both bathymetries). Cost on one T4 (gpu-dev).
Env of the base configuration: `matrix/calenv.sh`. Figure: `calibration_surface.png`.

## Step 0: minimum depth 5 m, no isolated-cell removal, no velocity clip

**Reproduction.** Native (ETOPO) depth, MINIMUM_DEPTH=5, ISOLATED_DEPTH_RATIO=0, Δt 10 min, χ 1.5 (old C field), SH:
it does **not** blow up (60 days, jobs 3948/3949), but the 5 m cell at (640, 50) (300.1ᵒE 62.3ᵒS) builds a steady
set-up of 13 m (η/H ≈ 2.6) with |u| 6.5 m/s, growing from day ~5 and saturating by day ~30. The same steady
mean flow (several metres of mean η, steady along-slope currents up to 12 m/s in 3000 m of water, ~2000 cells
> 0.5 m/s) exists in every earlier run near Antarctica (Ross/Amundsen, South Shetlands), also with 20 m minimum
depth and in the MOM6-depth configuration. With Cd = 0 a 5 m cell reaches η ≈ 25 m by day 5.

**Suspects, one at a time** (5–30 day runs):
- (b) bottom drag: `ocean_simulation` default is `implicit_bottom_drag = true` (IMEX flux, λ = −Cd|u|); barotropic_tides.jl
  did not override it. In the split-explicit solver the implicit step is applied to u before the slow tendency is
  formed as the realised (u − u⁻)/Δt, so the drag does reach the barotropic mode. New option `BOTTOM_DRAG_FORM=explicit`:
  identical results (max η 4.65 vs 4.54 m at day 5; 13.5 vs 13.3 m at day 31). Not the cause.
- (c) η/H: the large η/H is the symptom (mean set-up), not the cause.
- forcing switched off at day 12 (`M2_FORCING_STOP_DAY`): the mean set-up **keeps growing** to 12.3 m at day 30,
  identical to the forced run → a self-excited linear instability seeded by the spin-up, not tidal rectification.
- Coriolis off (`CORIOLIS=off`, diagnostic): mean η exactly 0.
- Coriolis schemes EnstrophyConserving / EnergyConserving / ActiveWeighted*: all show it.
- (a) thickness-weighted vs flat Laplacian, (d) depth steps: they matter only through the Coriolis term (below).

**Root cause: Oceananigans' Coriolis is not energy conserving over partial cells.** `EnergyConserving` is
Fᵘ = ℑʸ(f ℑˣ(A_y v)) / A_y^u (A_y = Δx Δz). Summed by parts, the work of the two components pairs
ℑˣ(Δx Δz v) ℑʸ(Δy u) against ℑʸ(Δy Δz u) ℑˣ(Δx v); these cancel only for uniform Δz. With true-depth partial cells
(5 m next to 4000 m) Coriolis does net work and a steady geostrophic mode grows on steep slopes where f is large,
until bottom drag balances it. The enstrophy and active-weighted variants share the mismatch. Proof:

| test (native depth, 5 m min, unforced, inviscid, no drag, 0.05 m/s white-noise seed, Δt 10 min) | E(day 0.25) | E(1) | E(2) | E(3) | E(4) | E(8) |
|---|---|---|---|---|---|---|
| EnergyConserving | 1.619e15 | 1.825e15 | 2.490e15 | 4.070e15 | 9.627e15 | NaN (day ~4.5) |
| ThicknessWeightedEnergy | 1.602e15 | 1.596e15 | 1.595e15 | 1.594e15 | 1.594e15 | 1.5935e15 |

(E = ½ΣV|u|² + ½gΣA_z η², m⁵ s⁻²; the first step removes the divergent part of the seed.) Algebraic check on a random
field (`coriolis_work.jl`): the work summed onto (Face, Face) corners is exactly 0 at every corner for the new form;
for EnergyConserving Σ|corner work| = 1.7e12, net −1.1e10, concentrated at partial-cell steps. (A per-face map is
not zero for any scheme: work is exchanged between neighbouring faces; only the corner sums cancel.)

**Fix: `CORIOLIS_SCHEME=ThicknessWeightedEnergy`** (forcing in barotropic_tides.jl), the Sadourny energy form on
transports with corner potential vorticity f/h:
Fᵘ = ℑʸ(f ℑˣ(A_y v)/h) A_x/V,  Fᵛ = −ℑˣ(f ℑʸ(A_x u)/h) A_y/V,
whose work cancels exactly per corner for any h > 0. h = mean of the wet face thicknesses Δzᶠᶜᶜ(i, j), Δzᶠᶜᶜ(i, j−1),
Δzᶜᶠᶜ(i, j), Δzᶜᶠᶜ(i−1, j) (default; `CORIOLIS_CORNER_THICKNESS=cells` uses the mean of the wet cells instead: also
energy-conserving, but biased against the min-of-neighbours face thickness on slopes, weakening Coriolis there:
EOT20 deep 7.58 vs 6.28 cm, phase offset −8.4ᵒ vs −4.0ᵒ). Reduces to EnergyConserving for uniform Δz; f at
(Face, Face); land handled by wet masks; the fold row through the halo fill.
Result: mean η ≤ 0.005 m everywhere, no cell with mean speed > 0.1 m/s (s0_twe); native 5 m / no removal 60 days stable
(job 3984, max|u| 1.2 m/s). No nonlinear thickness or wetting/drying needed. Worth an Oceananigans issue (the same
holds for any `HydrostaticSphericalCoriolis` scheme on `PartialCellBottom` grids).

Open: Lincoln Sea resonance, M2 amplitude 5–6 m at (723, 454) 295ᵒE 81ᵒN (real ~0.3 m), unchanged by the fix. It
sits in a **one-cell-wide channel** along row j = 454 (i ≈ 717–724, 61–376 m deep, rows 453 and 455 land), 18 rows
south of the fold: a grid-scale closed channel resonating, not a fold artefact.

## Step 1: biharmonic-only viscosity

`VISCOSITY=biharmonic`, `BIHARMONIC_TIMESCALE=τ` (days): A = Δ⁴/τ, Δ² = 2Δx²Δy²/(Δx²+Δy²), flat (not thickness-
weighted) 2D operator, free slip, no Laplacian, no Smagorinsky, bounded A ≤ Δ⁴/(64Δt) (the bound is hit
everywhere or nowhere: iff τ < 64Δt = 0.44 d at Δt 10 min, 0.22 d at 5 min; logged). Smoke tests at Δt 10 min,
MOM6 depth 5 m: τ = 1, 3, 10 d all stable over 5 days (max|η| 3.8 m, max|u| 1.9–2.0 m/s at Antarctic/Arctic hotspots);
τ = 0.5 d stable for 60 days at Δt 5 min. Old configurations are unchanged (BIHARMONIC_TIMESCALE alone still uses
the old closure).

## Time step and the split-explicit filter

No-SAL twins vs the differentiable session's explicit model (Δt 22.9 s), old Coriolis on both sides, deep:

| Δt | amplitude ours / explicit | rms|ΔZ|/√2 | mean |Δphase| |
|---|---|---|---|
| 10 min | 0.973 | 1.17 cm | 1.22ᵒ |
| 5 min | 0.986 | 0.65 cm | 0.77ᵒ |
| 2.5 min | 0.993 | 0.43 cm | 0.62ᵒ |

The amplitude deficit is linear in Δt: the time-averaging filter of the split-explicit free surface damps M2.
At 10 min, the cosine averaging kernel recovers about half of the 10 → 2.5 min change (as Δt 5 min would);
doubling the substep count (FREE_SURFACE_CFL 0.35) changes nothing; the constant kernel is wrong (17% error).
New options `FREE_SURFACE_KERNEL`, `FREE_SURFACE_CFL`. Chosen: **Δt 5 min, default kernel**; residual
~1.4% amplitude deficit (~0.65 cm deep) relative to the explicit solution, absorbed by χ. Check at the optimum:
Δt 3 min gives 6.15 vs 6.17 cm (deep EOT20) — negligible.

## Step 2: sweep (Δt 5 min, MOM6 aw depth, min 5 m, JSL_RH × χ semi-implicit, Cd 0.003, g 9.8, M2 only, ThicknessWeightedEnergy, SH n=40 PREM, 60 days, fit 30–60)

| χ | τ (d) | EOT20 deep (tuned) | amp / phase | offset ᵒ | TPXO9 deep (held out) | gauges all / > 200 m (held out) |
|---|---|---|---|---|---|---|
| 0.6 | 30 | 7.19 | 4.29 / 5.77 | −6.3 | 7.47 | 50.6 / 19.9 |
| 0.8 | 30 | 6.35 | 3.81 / 5.08 | −4.4 | 6.61 | 49.7 / 19.1 |
| 0.9 | 30 | 6.19 | 3.81 / 4.87 | −3.6 | 6.42 | 49.3 / 18.8 |
| 1.0 | 30 | 6.17 | 3.94 / 4.75 | −2.8 | 6.38 | 49.0 / 18.6 |
| 1.1 | 30 | 6.27 | 4.17 / 4.68 | −2.1 | 6.45 | 48.8 / 18.4 |
| 1.2 | 30 | 6.46 | 4.46 / 4.68 | −1.4 | 6.62 | 48.6 / 18.3 |
| 1.4 | 30 | 7.04 | 5.18 / 4.77 | −0.2 | 7.15 | 48.3 / 18.3 |
| 1.0 | 100 | 6.17 | 3.94 / 4.75 | −2.8 | 6.38 | 49.1 / 18.6 |
| 1.0 | 10 | 6.17 | 3.94 / 4.75 | −2.8 | 6.38 | 48.9 / 18.5 |
| 1.0 | 3 | 6.17 | 3.93 / 4.75 | −2.8 | 6.38 | 48.5 / 18.3 |
| 1.0 | 1 | 6.15 | 3.90 / 4.75 | −2.8 | 6.36 | 47.6 / 17.8 |
| **1.0** | **0.5** | **6.13** | 3.88 / 4.75 | −2.7 | **6.34** | **46.5 / 17.3** |
| 0.9 | 1 | 6.17 | 3.77 / 4.88 | −3.5 | 6.40 | 47.8 / 18.0 |
| 1.1 | 1 | 6.26 | 4.14 / 4.70 | −2.0 | 6.44 | 47.4 / 17.7 |

(14 runs. Earlier 10-min runs with the cell-thickness Coriolis, superseded: χ 0.6–1.2 → 8.26, 7.58, 7.37, 7.52.)

**No τ–χ ridge.** The deep error depends on χ only: a parabola with minimum at χ ≈ 1.0 (0.9–1.1 within 0.1 cm).
The viscosity is irrelevant for deep M2 over τ = 0.5–100 d (≤ 0.04 cm); it matters only for coastal gauges
(> 200 m: 18.6 → 17.3 cm as τ goes 100 → 0.5 d), where it damps grid-scale shelf/channel responses.

**Optimum: χ = 1.0 (JSL_RH unscaled), τ = 0.5 d** — deep EOT20 6.13 cm, held-out TPXO9 6.34 cm, gauges
46.5 / 17.3 cm. Physically: the JSL linear wave drag as computed (no tuning factor) is optimal once the spurious
Coriolis energy source is gone; the strongest stable biharmonic damping (τ ≈ 2× its stability bound) is preferred
only by the held-out gauges, so τ is weakly constrained (anything 0.5–3 d is equivalent in the deep ocean).
The deep error splits as 3.9 cm amplitude, 4.75 cm phase, mean phase offset −2.7ᵒ (model lags).

## Step 3: final checks at χ 1.0, τ 30 d

| run | EOT20 deep | offset | TPXO9 deep | gauges all / > 200 | wall s / model day |
|---|---|---|---|---|---|
| Δt 5 min (SAL at every substage, default) | 6.17 | −2.8ᵒ | 6.38 | 49.0 / 18.6 | 34.6 |
| Δt 5 min, SAL once per step (`SAL_UPDATE_INTERVAL=1`) | 6.22 | −2.4ᵒ | 6.41 | 48.9 / 18.5 | 25.9 |
| Δt 3 min | 6.15 | −2.9ᵒ | 6.37 | 49.1 / 18.6 | 52.2 |
| convolution PREM (Float32 accurate), every step, Δt 5 min | running (job 4029) | | | | ~1016 (≈ 17 h for 60 days) |

SH SAL update 22 ms; convolution 3.3 s per update (Float32 accurate; Float64 50–57 s, infeasible on T4).

## Changed files (not committed)

- `GlobalKiloMischief/tides/barotropic_tides.jl`: `CORIOLIS_SCHEME=ThicknessWeightedEnergy` (+ `CORIOLIS_CORNER_THICKNESS`),
  `VISCOSITY=biharmonic`, `BOTTOM_DRAG_FORM`, `FREE_SURFACE_KERNEL`, `FREE_SURFACE_CFL`, `SAL_UPDATE_INTERVAL=1` = once
  per step, diagnostics `CORIOLIS=off`, `M2_FORCING_STOP_DAY`, `SEED_VELOCITY`, `ENERGY_DIAGNOSTIC`. Defaults unchanged.
- `analysis/matrix/`: `calenv.sh`, `coriolis_work.jl`, `calib_figure.jl`, `calib_scores.tsv`, `calibration_surface.png`,
  `dtcompare.sh`, `fscompare.sh`, `CALIBRATION.md`.
