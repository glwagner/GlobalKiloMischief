# GlobalKiloMischief

Global, kilometer-scale ocean–sea ice simulations with [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl)
and [NumericalEarth.jl](https://github.com/NumericalEarth/NumericalEarth.jl), aiming to reproduce the MITgcm LLC4320
simulation (1/48°, 90 levels, tides) on a `TripolarGrid`. We start at 1/6° on one GPU and 1/12° on four GPUs, then
refine. Runs on NCSA DeltaAI (GH200).

| Script | Configuration |
|---|---|
| `smoke_test.jl` | Full 1/6° configuration, 1 GPU, 300 steps. Downloads and caches all inputs. |
| `sixth_degree.jl` | 1/6° (2160 × 1020 × 100), 1 GPU |
| `twelfth_degree.jl` | 1/12° (4320 × 2040 × 100), 4 GPUs on one node (`NCCLDistributed`, `Partition(1, 4)`) |
| `download_data.jl` | Fetches the forcing for a window on the login node (CPU only) |
| `setup.jl` | Shared grid, model, and simulation constructors |

## Configuration

- **Grid**: `TripolarGrid` from 80°S, 100 z⋆ levels (1.4 m at the surface, 320 m at 6000 m), ETOPO2022 bathymetry
  (`GridFittedBottom`, active cells map). The bathymetry is built on the CPU for the whole grid and partitioned across
  ranks. Only the connected world ocean is kept, so the Black and Caspian Seas are land. `connect_basins!` keeps
  the Mediterranean: if the Gulf of Cádiz and the Alboran Sea fall in different basins, it deepens to the 280 m
  Camarinal Sill the cells along the least-deepening, face-connected path between them. In practice the strait is
  closed at 1°, 1/3° and 1/6° (2–4 cells deepened) and already open at 1/2°, 1/4° and 1/12°.
- **Ocean**: `ocean_simulation` defaults: `WENOVectorInvariant` momentum, `WENO(order=7)` tracers (both
  vertically implicit where the vertical CFL exceeds 0.5), `CATKEVerticalDiffusivity`, `SplitRungeKutta3`,
  `SplitExplicitFreeSurface` with a substep count sized for `max_Δt`. No GM, no horizontal viscosity. One change:
  near walls the tracer WENO falls back to first-order upwind (`minimum_buffer_upwind_order = 1`) instead of
  second-order centered. Centered fluxes leave bottom corner cells undamped, and overflow cells in the Faroe Bank
  Channel and Denmark Strait cooled to -12.7 ᵒC in 10 days.
- **Sea ice**: ClimaSeaIce, zero-layer thermodynamics with snow, mEVP with 120 substeps (relaxation parameter ≤ 120),
  incremental remapping of ice volume and concentration (`IncrementalRemapping()`, forward Euler), which keeps thick ice
  from piling into coastal cells at low concentration, as `WENO` advection of the two separately did.
- **Initial condition**: GLORYS12 daily T, S, sea ice thickness and concentration on 2015-01-01. Ocean at rest.
- **Forcing**: ERA5 hourly single levels (10 m wind, 2 m temperature and dewpoint, surface pressure, precipitation,
  downwelling short and longwave), JRA55-do rivers and icebergs. Latitude-dependent ocean albedo.
- **Time step**: adaptive, for maximum throughput. A `TimeStepWizard` limits the horizontal advective CFL to 0.5
  (vertical advection is implicit) and grows Δt by at most 5% every 10 steps, from 1 min (1/6°) or 30 s (1/12°) up to
  `max_Δt`, a script argument (default 20 min at 1/6°, 10 min at 1/12°). The barotropic substep count is sized for
  `max_Δt`, so every Δt ≤ `max_Δt` is stable for the free surface. Advection permits about an hour at 1/6°, so the
  practical limits are internal gravity waves on the smallest cells (near Antarctica) and the coupled sea ice; find
  them by raising `max_Δt` until a run fails. Each `max_Δt` writes to its own run directory.
- **Throughput**: every 50 steps the log reports Δt, wall time per step, and SYPD (simulated years per wall-clock
  day) over those steps, together with extrema and GPU memory.
- **Output** (`$RUN_DIRECTORY/<name>/`, on `/work/nvme`): daily surface T, S, e, u, v, w; η; sea ice h, ℵ, u, v. Checkpoints every
  5 simulated days, so 1-hour interactive jobs of 5 days each can be chained. Each run picks up from the latest checkpoint, so you extend a run by resubmitting with a
  larger `stop_days`.

## Pinned versions

Julia ≥ 1.12.3 is required: earlier versions assume 4 KiB pages on aarch64, and JLD2's memory-mapped writes then
fail on DeltaAI's 64 KiB pages (`SystemError("msync", 22)`, JLD2 #702).

`Project.toml` pins NumericalEarth at `9ef864a` (branch `glw/global-kilo-mischief`), which is `main` plus four open PRs:

- #742: fixes GLORYS inpainting. Deep levels and land were left at zero, and atoll columns got surface water copied down to the seafloor.
- #743: lets ERA5 supply mean sea level pressure instead of surface pressure.
- #769: builds JRA55 river routing on distributed grids (it crashed rebuilding the grid on the CPU).
- #770: adds `frazil_formation_depth`. We set 50 m: deep cells cooled below freezing by advection undershoots at steep
  topography (900 m near the Faroes) made up to 5.5 m/day of surface ice over 9 ᵒC water.

Oceananigans is pinned at `66e3dda` (branch `glw/v0.113.6-north-fold`): the v0.113.6 release plus #6208, which
drops duplicated north-fold work for y-partitioned distributed tripolar grids, as in the 4-GPU `Partition(1, 4)` run.
ClimaSeaIce is 0.5.11, which already includes the distributed and tripolar EVP fixes (#137, #152, #146, #135) and the
metric-aware stresses (#164).

Oceananigans 0.114 also reworks the active-cell map (#6180), which ClimaSeaIce #184 builds on. It needs ClimaSeaIce and
NumericalEarth to accept 0.114 first.

Upstream PRs to watch: Oceananigans #6111 (SSPRK3), #6082 (CUDA graphs for the barotropic solver), #6152
(implicit drag with split RK) and #5888 (active-cell load balancing); ClimaSeaIce #184 (sea ice on its own
immersed grid) and #180 (CUDA graphs for EVP).

## Running

```bash
# One-time: export credentials in ~/.bashrc (slurm/env.sh sources it)
#   COPERNICUSMARINE_SERVICE_USERNAME, COPERNICUSMARINE_SERVICE_PASSWORD     GLORYS
#   CDSAPI_URL=https://cds.climate.copernicus.eu/api, CDSAPI_KEY=<token>    ERA5 (accept the single-levels licence)
sbatch slurm/precompile.sbatch

source slurm/env.sh && julia --project download_data.jl 60   # login node: bathymetry + forcing for 60 days
sbatch slurm/smoke_test.sbatch
sbatch slurm/sixth_degree.sbatch 10 20                       # stop_days max_Δt_minutes; rerun with more days to extend
sbatch slurm/twelfth_degree.sbatch 10 10
```

## Known limitations

- JRA55-do runoff ends on 2019-12-31, which caps hindcasts started in 2015 at five years.
- No tides yet. They are planned once the runs are stable, as an equilibrium tidal potential through the
  pluggable barotropic potential of NumericalEarth #736.
