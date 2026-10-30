# Tides at 1ᵒ

Global tides forced by the astronomical tidal potential, starting with a single-level (shallow-water) ocean on the
1ᵒ tripolar grid. This directory has its own environment:

```
julia --project=tides -e 'using Pkg; Pkg.instantiate()'
julia --project=tides tides/single_level.jl <CPU | CUDA | Metal> <Float64 | Float32> <Δt in minutes> <days>
```

For Metal, `Pkg.add("Metal")` in this environment first; it is not a dependency so that the environment installs on Linux.

## What each package provides

| Package | Version | Provides |
|---|---|---|
| [Tidejinks](https://github.com/glwagner/Tidejinks.jl) | main 2539717 | `TidalPotential(grid, dates; reference_frame)`: the degree-two potential of the Sun and Moon, with the Love number factor 1 + k₂ − h₂ = 0.69, from SPICE ephemerides. It matches MITgcm's `tides_calc.F` to 2e-11. It is evaluated in kernels from five coefficients sampled hourly (error 4e-9 in Float64, 1.5e-6 in Float32). |
| [NumericalEarth](https://github.com/NumericalEarth/NumericalEarth.jl) | branch `glw/gkm-tides-stack`, cc4230cc | `ocean_simulation(grid; barotropic_potential)` from NumericalEarth#736, which forces `u` and `v` with `-∇Φ` evaluated at each Runge-Kutta stage time. The branch is main plus the open PRs #736, #769, #770, #773, #774, #779 and #800. |
| [ClimaSeaIce](https://github.com/CliMA/ClimaSeaIce.jl) | main 38a99e5 (0.6.0) | Required by NumericalEarth. 0.6 is needed for distributed runs with Oceananigans 0.114. |
| [Oceananigans](https://github.com/CliMA/Oceananigans.jl) | 0.114.1 | The hydrostatic model, z⋆ coordinate, split-explicit free surface. |
| SPICE, JLD2, CairoMakie, CUDA | registered | Ephemerides, output, plots, NVIDIA GPUs. |

Tidejinks downloads the NAIF kernels (DE440, leap seconds, Earth frames, GM) to a scratch space when it loads, and
NumericalEarth downloads ETOPO2022 bathymetry. Compute nodes without internet need both cached from a login node first,
by running `single_level.jl` there for a fraction of a day on the CPU.

## The single-level configuration

`single_level.jl` builds:

- a 360 × 170 tripolar grid with one z⋆ level and ETOPO2022 bathymetry (minimum depth 10 m, minor basins removed),
  with a partial-cell bottom whose minimum fractional height is 1e-3, so that every column has its true depth;
- `ocean_simulation` with uniform temperature and salinity, no closure, no radiation, implicit quadratic bottom drag
  (Cᴰ = 0.003), and the tidal potential in the `IAU_EARTH` reference frame (as in MITgcm) from 2015-01-01;
- hourly output of η to `tides/output`.

Environment variables select variants: `MOMENTUM_ADVECTION=none | centered` (default `WENOVectorInvariant()`),
`HORIZONTAL_VISCOSITY` (m² s⁻¹), `NO_SLIP=true`, `IMPLICIT_BOTTOM_DRAG=false`, `MINIMUM_DEPTH`, `RUN_SUFFIX`.

## Internal tides

`internal_tides.jl <CPU | CUDA | Metal> <Float64 | Float32> <cells per degree> <Nz> <Δt in minutes> <days>` follows
Simmons, Hallberg & Arbic (2004): an ocean at rest with horizontally uniform stratification N = N₀ exp(z / d)
(N₀ = 5.2e-3 s⁻¹, d = 1300 m) forced only by the tidal potential. It saves hourly the isopycnal displacement ζ at about
1000 m, the vertical velocity at about 200 m, and η. A biharmonic viscosity ν = Az² / τ (`BIHARMONIC_TIMESCALE`, default
30 days) removes the wave energy that cascades to the grid scale; without it the field fills with grid-scale noise
within a week. At 1/2ᵒ the mode-1 M2 wavelength (about 100 km) is only two
cells; use 1/4ᵒ or finer.

Analysis: `harmonic_analysis.jl` fits M2, S2, K1, O1 over days 10–30 and compares runs; `animate_globe.jl <run> <η | ζ | w> <first day> <last day>` animates a
field on the sphere (`HIGH_PASS=true` subtracts the 25-hour running mean); `locate_blowup.jl` finds where η grows; `spin_down.jl` checks that bottom drag slows the depth-mean flow.

## Results so far (Metal, CPU; 2026-10-22)

- With momentum advection, every run blows up between day 2 and day 4: shallow coastal cells drain until the column
  dries. The blow-up comes sooner as Δt decreases (4.2 days at 20 min, 2.15 days at 1 min), so it is not a time-stepping
  instability. Viscosity up to 5e4 m² s⁻¹, no-slip walls, and a 50 m minimum depth do not prevent it.
- Without momentum advection (`MOMENTUM_ADVECTION=none`) runs are bounded but resonate to 17 m on the shelf off the
  Antarctic Peninsula. The global rms M2 amplitude is 0.66 m, about twice the observed value, as expected without
  self-attraction and loading or internal-wave drag.
- Metal Float32, CPU Float32 and CPU Float64 agree: at equal Δt the M2 difference is 7e-4 of the rms amplitude.
- M2 converges at first order in Δt (relative differences 8e-2, 3e-2, 2e-2 at Δt = 20, 10, 2.5 min against 5 min).
