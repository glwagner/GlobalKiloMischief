# GlobalOceanSeaIce

Eddying global ocean–sea ice simulations with [NumericalEarth.jl](https://github.com/NumericalEarth/NumericalEarth.jl)
on a `TripolarGrid`. Runs on NCSA DeltaAI (GH200).

| Script | Configuration |
|---|---|
| `smoke_test.jl` | Full 1/6° configuration, 1 GPU, 300 steps. Downloads and caches all inputs. |
| `sixth_degree.jl` | 1/6° (2160 × 1020 × 100), 1 GPU |
| `twelfth_degree.jl` | 1/12° (4320 × 2040 × 100), 4 GPUs on one node (`NCCLDistributed`, `Partition(1, 4)`) |
| `download_data.jl` | Fetches the forcing for a window on the login node (CPU only) |
| `setup.jl` | Shared grid, model, and simulation constructors |

## Configuration

- **Grid**: `TripolarGrid` from 80°S, 100 z⋆ levels (1.4 m at the surface, 320 m at 6000 m), ETOPO2022 bathymetry
  (`GridFittedBottom`, active cells map).
- **Ocean**: `ocean_simulation` defaults only: `WENOVectorInvariant` momentum, `WENO(order=7)` tracers (both
  vertically implicit where the vertical CFL exceeds 0.5), `CATKEVerticalDiffusivity`, `SplitRungeKutta3`,
  `SplitExplicitFreeSurface` with a substep count sized for `max_Δt`. No GM, no horizontal viscosity.
- **Sea ice**: ClimaSeaIce, zero-layer thermodynamics with snow, mEVP with 120 substeps (relaxation parameter ≤ 120),
  `WENO(order=7)` advection.
- **Initial condition**: GLORYS12 daily T, S, sea ice thickness and concentration on 2015-01-01. Ocean at rest.
- **Forcing**: ERA5 hourly single levels (10 m wind, 2 m temperature and dewpoint, surface pressure, precipitation,
  downwelling short and longwave), JRA55-do rivers and icebergs. Latitude-dependent ocean albedo.
- **Time step**: adaptive. `TimeStepWizard` with horizontal CFL 0.5, starting from 1 min (1/6°) or 30 s (1/12°)
  and growing by at most 5% every 10 steps up to `max_Δt` (10 min at 1/6°, 5 min at 1/12°). Vertical advection
  is implicit, so it is left out of the CFL.
- **Output** (`$RUN_DIRECTORY/<name>/`): daily surface T, S, e, u, v, w; η; sea ice h, ℵ, u, v. Checkpoints every
  10 simulated days. Each run picks up from the latest checkpoint, so you extend a run by resubmitting with a
  larger `stop_days`.

## Pinned versions

`Project.toml` pins NumericalEarth `d07eb24` and Oceananigans `main` at `ac83d85` (still 0.113.5). The `main` pin
picks up changes that have not been released yet:

- #6154: the serial tripolar fold no longer needs extended halos.
- #5897 and #6081: fewer host–device syncs and allocations in distributed halo exchange.

ClimaSeaIce 0.5.11 already includes the distributed and tripolar EVP fixes (#137, #152, #146, #135) and the
metric-aware stresses (#164).

Upstream PRs to watch: Oceananigans #6111 (SSPRK3), #6082 (CUDA graphs for the barotropic solver), #6152
(implicit drag with split RK) and #5888 (active-cell load balancing); ClimaSeaIce #184 (sea ice on its own
immersed grid) and #180 (CUDA graphs for EVP).

## Running

```bash
# One-time: credentials
#   ~/.cdsapirc                   ERA5 (CDS personal access token; accept the ERA5 single-levels licence)
#   ~/.copernicusmarine.env       export COPERNICUSMARINE_SERVICE_USERNAME=... / _PASSWORD=...
sbatch slurm/precompile.sbatch

source slurm/env.sh && julia --project download_data.jl 60   # login node: forcing for 60 days
sbatch slurm/smoke_test.sbatch
sbatch slurm/sixth_degree.sbatch 10                          # then 30, 90, ... (each picks up)
sbatch slurm/twelfth_degree.sbatch 10
```

## Known limitations

- `ERA5PrescribedAtmosphere` keeps the derived specific humidity in memory for the whole forcing window
  (about 100 MB per simulated day on each GPU). Runs longer than a few months need NumericalEarth to derive
  humidity from dewpoint at load time, so that it is windowed like the other ERA5 fields.
- JRA55-do runoff ends on 2019-12-31, which caps hindcasts started in 2015 at five years.
