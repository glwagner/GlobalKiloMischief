# MOM6-matched and native configurations

Below are the two configurations as the code actually runs them. I took them from `barotropic_tides.jl`, the two config scripts (`matrix/m6env.sh` for MOM6-matched, `matrix/bestenv.sh` for native), `RESULTS.md`, and the NumericalEarth `ocean_simulation` defaults. Anything I couldn't confirm in code is marked **[unverified]**.

## 1. Governing equations (shared form)

Both configurations run the same single-layer, linear shallow-water system on a 0.36° tripolar grid (1000 × 472 cells). There is one layer of fixed thickness, depth H(x) is set by `PartialCellBottom` with minimum fractional cell height ε ≤ 10⁻³ (so true depth is kept), and there is no momentum advection (`MOMENTUM_ADVECTION=none`).

**Momentum:**

∂ₜ**u** + f ẑ×**u** = −g∇η − ∇Φ_tide − ∇Ψ_SAL + **D**_visc − (C_d |**u**|/H) **u** − r(x) **u**

**Continuity:**

∂ₜη + ∇·(H**u**) = 0

- **Depth in continuity:** the code uses fixed H. MOM6 uses D + η in the continuity equation and in the drag, so it is nonlinear there. Ours is linear in η **[unverified: I believe Oceananigans' split-explicit transport uses H, not H + η, here]**.
- **Bottom drag:** C_d = 0.003, quadratic, applied as a bottom flux on the single cell, so the acceleration is C_d|u|u/Δz with Δz = H. The `ub` background velocity inside |u| is set by NumericalEarth **[unverified value]**.
- **Wave drag:** r = χ C(x)/H, the Jayne–St. Laurent linear form, averaged onto u/v faces. When `WAVE_DRAG_IMPLICIT=true` it is written as r/(1 + rΔt); otherwise it is capped at 1/(2Δt).
- **SAL:** Ψ = −Σₙ γₙ pᵇₙ with γₙ = (3ρ₀/ρₑ)(1 + k′ₙ − h′ₙ)/(2n+1) and PREM Love numbers (`LOADING_LOVE_NUMBERS=prem`). The matrix runs this as SH with n=40, or as convolution, or replaces it with scalar Ψ = −0.1 g η.
- **Time stepping:** split-explicit free surface with CFL 0.7 substeps and `SplitRungeKutta3` baroclinic stepping. That is NumericalEarth's default; neither script overrides it.

## 2. Side by side

| | MOM6-matched (`m6env.sh`) | Native (`bestenv.sh`) |
|---|---|---|
| Bathymetry | Chen et al. p04arc9 `aw` depth, nearest-cell binned onto our grid (`regrid_mom6.jl`). Max 9314 m, so the layer is as deep as the deepest cell. | NumericalEarth `regrid_bathymetry` (ETOPO), 3 interpolation passes, major basins only, **clipped at 6000 m by the script's own `depth = 6000` (since fixed: regrid on an 11000 m grid, layer as deep as the deepest cell)** |
| Minimum depth | 5 m. 744 cells MOM6 keeps wet above sea level become 5 m. | 20 m |
| Isolated cells | Kept (`ISOLATED_DEPTH_RATIO=0`) | Shallow cells with every wet neighbour ≥5× deeper become land, repeated until none are left |
| Wave drag field C | MOM6 `JSL_RH`, regridded | Same `JSL_RH` |
| Wave drag scale χ | 0.8 (Table 1's SH value) | 1.1 |
| Wave drag H | MOM6 depth, consistent with how `JSL_RH` was derived | **Our** depth, a different H from the one `JSL_RH` was built on |
| Wave drag stepping | Semi-implicit. Max r = 1.02×10⁻⁵ s⁻¹, rΔt = 0.0018. | Explicit with cap, though nothing reaches the cap |
| Viscosity | MOM6 form as a flat 2D forcing. Laplacian K_H = min(0.03 Δ, Δ²/8Δt), at most 1234 m²/s. Biharmonic A_H = max(0.01 Δ³, 0.06 Δ⁴\|S\|), bounded by Δ⁴/64Δt. Δ² = 2Δx²Δy²/(Δx²+Δy²), free slip at land. | Uniform Laplacian ν = 1000 m²/s through Oceananigans' `HorizontalScalarDiffusivity`, which is thickness-weighted. No biharmonic. |
| Gravity g | 9.8 | 9.80665 (default) |
| Time step Δt | 3 min (ARGS) **[unverified: from RESULTS.md; I didn't see the submit line]** | 10 min |
| Forcing Φ_tide | M2 only: −0.69 × 0.2423 m × g cos²φ cos(ωt + 2λ + phase), nodal f and u fixed at 2015-01-16. The 0.69 is 1 + k₂ − h₂ for body tides. | Full Tidejinks potential, all constituents. Whether it includes the 1 + k − h factor is **[unverified]**. |
| Coriolis | `EnergyConserving` (MOM6 uses Sadourny 1975 energy) | `EnstrophyConserving` (default) |
| Bottom drag | C_d = 0.003. In MOM6, HBBL = 10 m only sets the sampled velocity in one layer. | C_d = 0.003 |
| Not matched to MOM6 | MAXVEL 6 m/s clipping, BOUND_CORIOLIS, BT_NONLIN_STRESS, split-RK2 stepping, D + η thickness | – |
| Run and fit | 60 days, fit over days 30–60 | Same |

## 3. Why native is set up the way it is

| Choice | Motivation | Status |
|---|---|---|
| ETOPO via `regrid_bathymetry` | It's the standard NumericalEarth path, and the depth our other models use. | The 6000 m clip came from the script's `depth = 6000`, not from NumericalEarth. |
| 20 m minimum, isolated-cell removal | Stability. Isolated 20 m cells in deep water went unstable, for example at 172°E 6°N. | Chosen for stability, not physics. |
| No advection, Laplacian 1000 m²/s, no biharmonic | WENO advection blew up next to 20 m cells. Oceananigans' biharmonic blows up at partial-cell steps because it's thickness-weighted. | Sweep showed ν doesn't matter (0.01 cm). |
| `JSL_RH` × 1.1 | The physically computed field (from 1/12° topography) at scale about 1. It beat our 0.5° Egbert–Ray field and ranked the same on held-out data. | Tuned on deep EOT20, but only on our depth. |
| Δt 10 min | Cheap. Δt of 5 and 2.5 min changed the score by 0.1 cm or less. | Tested |
| Full potential | Realism | M2-only tested: no change |

## 4. What explains native's larger deep error (8.31 vs 6.04 cm)

Ruled out by earlier tests in the native setup:
- **Δt:** 2.5–10 min gave 8.56–8.69.
- **M2-only vs full forcing:** 8.63 either way.
- **Viscosity form and size:** within 0.06 cm.
- **Coriolis scheme:** about 7×10⁻⁴ of the force.
- **g:** 0.07% in g is 0.03% in wave speed.
- **Implicit vs explicit drag:** 8.55 vs 8.56.

What's left is **the geometry and everything derived from it**, which I rank most likely:

1. **Bathymetry source and smoothing.** Tide speed is c = √(gH), and basin-scale M2 phase is the integral of 1/c along the wave's path. The native error maps show large-scale patterns in the central Pacific and South Atlantic, which is what a wave-speed or basin-geometry error produces. A dissipation error would show up locally instead. The extra interpolation passes and the 6000 m clip both change H. The clip affects little area, so I'd expect the source and smoothing to matter more **[inference]**.
2. **Shelf and coastline geometry:** the 20 m minimum, removed isolated cells and major-basins-only mask. Shelf resonances such as Hudson Strait, Patagonia and the North-West European shelf feed back on the open-ocean tide through energy drain and the near-resonant Atlantic response. This would also explain why native is better at gauges but worse in the deep ocean.
3. **Wave drag mismatch.** `JSL_RH` was computed for MOM6's topography, but native divides it by our H. Where our H is shallower or deeper, the drag rate is wrong. Native also needs χ = 1.1 where matched uses 0.8, which suggests compensation.

The test is to run native and change one thing at a time, using SH with n=40 (see section 6, runs 1–4).

## 5. How to improve on MOM6-matched

This setup is already at 6.04 cm, close to MOM6 itself at 6.35–6.58 cm. The improvements should be physical, not tuned to EOT20:

1. **Retune wave drag for this setup.** χ = 0.8 is MOM6's value, never optimised here. Sweep χ = 0.6, 0.7, 1.0 and 1.2, choose on deep EOT20, and validate on TPXO9 and gauges. This is the parameter with the strongest leverage.
2. **Use a wave drag field built from the actual model depth and stratification.** Recompute C = κh²N_b√(1 − f²/ω²)/2 using our H and an N_b climatology, so the drag field and the depth agree.
3. **Raise resolution.** Run 0.25° and 1/6° with the same inputs. At 0.25° native SH already improved from 8.63 to 8.19. Coastline and shelf representation is the biggest remaining approximation.
4. **Fix the open defects.** The tide-driven steady flow near Antarctica reaches 12 m/s, and the 5–6 m Arctic resonance sits at the fold. Both contaminate high latitudes and possibly global energy. Candidates are MOM6's velocity clip, a grid-scale mode at the fold, and partial-cell pressure gradients near steep slopes. These need a root cause, not a clip.
5. **Use D + η thickness in continuity and drag**, if Oceananigans doesn't already. On shelves where η/H is not small, it changes M2 phase and amplitude.
6. **Replace M2-only with the full potential.** It won't change M2 (tested), but it's needed for a realistic multi-constituent run and costs nothing.
7. **Keep SH n=40 as the default SAL.** Its accuracy equals convolution's at a fraction of the cost. If coastal detail matters, use a hybrid: SH plus a local convolution correction.

## 6. Proposed runs

All runs are 0.36°, 60 days, SH n=40 PREM. Each takes about 15–60 min on a T4.

| # | Run | Purpose |
|---|---|---|
| 1 | Native with MOM6 depth (`EXTERNAL_DEPTH=mom6`), native everything else | Isolates the bathymetry source |
| 2 | Native with 5 m minimum depth and `ISOLATED_DEPTH_RATIO=0` | Isolates shelf geometry |
| 3 | Run 1 with χ = 0.8 | Isolates the χ × H consistency of wave drag |
| 4 | Native with the 6000 m clip removed | Deep-ocean wave speed |
| 5–8 | MOM6-matched with χ = 0.6, 0.7, 1.0 and 1.2 | Tune the strongest parameter in the good setup |
| 9 | MOM6-matched at 0.25° | Resolution convergence |
| 10 | MOM6-matched with a self-consistent C on its own H | Physically derived drag |
| 11 | MOM6-matched with full potential and best χ | Production candidate |

These are 11 runs, around 6–8 GPU-hours, plus a few diagnostic runs for the Antarctic and Arctic defects.

## Unverified items

- Whether continuity uses H or H + η.
- The background velocity `ub` inside the quadratic drag.
- Whether Tidejinks applies the 1 + k − h body-tide factor in the native forcing.
- Δt = 3 min for MOM6-matched, which I took from `RESULTS.md` and didn't confirm from the submit line.
- That basin-scale phase error from bathymetry is the dominant cause, which is inference from error-map patterns; runs 1–4 test it.
