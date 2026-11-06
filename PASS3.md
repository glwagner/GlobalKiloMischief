# Pass 3: 1/12° global ocean–sea ice spin-up

90-day spin-up from GLORYS on 2021-06-30, forced by ERA5 hourly over 2021-06-30 to 2021-09-29, on 4 GPUs of one node.
Its first checkpoint and pickup double as the restart-at-CFL-1 gate test.

## Code and environment

- Branch `glw/pass3-ocean` of `github.com/glwagner/GlobalKiloMischief` (private); `Project.toml`, `Manifest.toml` and
  `LocalPreferences.toml` at the root. Julia 1.12.7.
- Oceananigans 0.116.1, ClimaSeaIce 0.7.1, NumericalEarth 0.9.1 at `4ae6a40005668518516cb22a96b134e0c44c3a92`.
  (The paper targets Oceananigans 0.117.0 once NumericalEarth and ClimaSeaIce release compat for it.)
- `setup.jl` (grid, bathymetry, model, output), `twelfth_degree.jl` (driver), `download_data.jl` (inputs).
- `slurm/*.sbatch` and `slurm/env.sh` are written for NCSA DeltaAI; on another cluster set the variables below
  instead of sourcing `env.sh`.

## Configuration

- `TripolarGrid` 4320 × 2040 × `Nz`, 80°S to the north poles, depth 6000 m, `ExponentialDiscretization` z⋆ levels
  (scale 1100 m). `Nz = 70` is the pass-3 option; 100 is the previous production setting.
- ETOPO2022 bathymetry, minimum depth 10 m, pits and dead ends filled, Gibraltar sill 280 m.
- `NCCLDistributed(GPU(); partition = Partition(1, 4))`: 4 ranks of 510 rows; NCCL halos, MPI scalar reductions.
- Adaptive Δt: `TimeStepWizard` on the horizontal advective CFL; the cap sizes the split-explicit substeps and must
  not change across a restart (a changed cap at pickup caused the CFL 1.0 blow-up in DeltaAI job 3498163).
  Δt is saved to `<run>/Δt` next to each checkpoint and restored on pickup.
- Checkpoints every 5 days (`cleanup = true`, about 65 GB each at Nz = 70), surface output every 6 h, daily η.

## Running

```
export START_DATE=2021-06-30
export RUN_DIRECTORY=/shared/greg_aeolus_earth/pass3_ocean/runs
export NUMERICALEARTH_DATA_DIRECTORY=<data root>   # see inputs below
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --project=. download_data.jl 90              # only fetches what is missing; also caches the bathymetry
export FI_PROVIDER=shm
mpiexec -n 4 julia --project=. twelfth_degree.jl 90 30 0.9 70   # stop day, Δt cap (min), CFL, Nz
```

Rerunning the same command picks up from the latest checkpoint. Keep the cap (30 min) the same for every
segment. The gate test is the first segment through day 5 and a pickup that holds CFL ≥ 0.8.
At CFL 1.0 the 1/6° run (DeltaAI job 3530121) reached Δt ≈ 60 min, so at 1/12° Δt should sit near 25–30 min
with the 30-min cap at or just above the CFL limit.

## Inputs (under `$NUMERICALEARTH_DATA_DIRECTORY`)

| Dataset | Variables | Window | Notes |
|---|---|---|---|
| GLORYS daily (`GLORYSDaily`, CMEMS GLOBAL_MULTIYEAR_PHY_001_030) | temperature, salinity, sea ice thickness, sea ice concentration | 2021-06-30 only | the last day of the reanalysis period; needs Copernicus Marine credentials to download |
| ERA5 hourly single levels (`ERA5HourlySingleLevel`) | 10 m u and v, 2 m temperature, 2 m dewpoint, mean sea level pressure, total precipitation, downwelling shortwave and longwave | 2021-06-30 to 2021-09-29 (stop + 1 day) | reported on falcon-dev `/data/coupled`; the file layout must match NumericalEarth's ERA5 naming or be re-fetched with CDS credentials |
| JRA55-do (`MultiYearJRA55`) | river and iceberg freshwater flux | 2021-06-30 to 2021-09-29 | check the JRA55-do version on hand covers 2021 |
| ETOPO2022 | bottom height | static | regridded once to the 4320 × 2040 grid and cached in `field_cache` |
