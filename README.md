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
| `download_data.jl` | Fetches the forcing for a window (CPU only; compute nodes reach CDS) |
| `setup.jl` | Shared grid, model, and simulation constructors |

## Results so far (2026-10-17)

| Config | Hardware | Δt | Step time | Throughput |
|---|---|---|---|---|
| 1/6° | 1 GH200 | ≤ 20 min | 1.57 s | 1.95 SYPD |
| 1/6° | 1 GH200 | ≤ 45 min, CFL 0.7 (avg ~37 min) | 1.7–1.8 s | **~3.4 SYPD** |
| 1/12° | 4 GH200 (1 node, NCCL) | 10 min (set limit) | 2.20 s (±1 ms) | 0.74 SYPD |
| 1/12° | 4 GH200 | ~20 min (projected) | ~2.2–2.3 s | ~1.3–1.5 SYPD |

- **1/6° stability:** 145 simulated days, stable. That run is now restarting with the corrected bathymetry.
- **1/12° status:** Δt above 10 minutes hasn't been tested yet, so the 20-minute row is a projection.
- **Weak scaling, 1 → 4 GPUs (on one node):** 71% efficiency per step, comparing 1/6° on 1 GPU with 1/12° on 4 GPUs.
- **MPI vs NCCL (1/12°, 4 GPUs, Δt = 10 min, 150 steps):** NCCL 2.19 s/step (±1 ms); CUDA-aware Cray MPICH 2.53 s/step
  (2.34–2.82 s), so NCCL is 13% faster and steadier. Cray MPICH works only with CUDA.jl's memory pool turned off
  (`JULIA_CUDA_MEMORY_POOL=none`: GTL's CUDA IPC cannot export pool memory); the pool makes no difference to NCCL.
- **Startup:** about 27 minutes per job, mostly slow reads from `/work/hdd`.
- **Weak scaling, 1 → 2 nodes (2026-10-18):** 1/12° (4320 × 2040) on 4 GPUs against 6120 × 2880 (≈ 1/17°) on 8 GPUs,
  same cells per GPU, Δt scaled with grid spacing (10 → 7 min), 160 steps:

  | Backend | 1 node, 4 GPUs | 2 nodes, 8 GPUs | Efficiency |
  |---|---|---|---|
  | NCCL (libfabric over Slingshot `cxi`, GPU Direct RDMA) | 2.19 s | 2.70 s (±3 ms) | **81%** |
  | Cray MPICH (CUDA-aware, memory pool off) | 2.53 s | 2.99 s (2.80–3.32) | 85% |

  1 GPU → 8 GPUs (1/6° → 1/17°) is 1.57 → 2.70 s/step, 58% overall. GPU memory reaches 82–91 GiB of 96 at 2.2M columns per GPU.

### Estimate: 1/48° on 64 GH200 (16 nodes)

- **Grid:** 17280 × 8160 × 100, about 14 billion cells. That's 2.2M columns per GPU, the same load per GPU as 1/6° on one GPU.
- **Step time:** 1/6° at matching substep counts runs about 1.75 s per step. At an expected 90% scaling efficiency,
  1/48° takes ~1.95 s per step.
- **Δt:** the advective CFL limit scales with grid spacing, so 37 min / 8 gives about 4.6 min (276 s).
- **Throughput:** 276 s / 1.95 s / 365 ≈ 0.39 SYPD. Reloads, output and checkpoints cost about 5–10% more, so
  **~0.35–0.37 SYPD**, about 4.3 simulated months per day.

| Assumption | Step time | SYPD |
|---|---|---|
| 90% efficiency, Δt 3.5–5 min | 1.95 s | 0.29–0.42 |
| Measured 8-GPU step time (2.70 s, 58% from 1 GPU), Δt 4.6 min | 2.7–3.0 s | **~0.25–0.28** |

The measured row assumes step time stops growing past 2 nodes: with a y partition each GPU exchanges halos with
two neighbors however many GPUs there are. An 8 × 8 partition has four neighbors per GPU, so it may cost more.

Caveats:

- **Sea-ice halos:** the sea-ice solver needs 243-row halos. A 64-way split in y alone leaves too few rows per GPU, so it
  needs an 8 × 8 split in x and y (or less wide halos). It hasn't been verified that `TripolarGrid` supports a
  two-dimensional split.
- **Unmeasured:** scaling beyond 2 nodes.

## Configuration

- **Grid**: `TripolarGrid` from 80°S, 100 z⋆ levels (1.4 m at the surface, 320 m at 6000 m), ETOPO2022 bathymetry
  (`GridFittedBottom`, active cells map). The bathymetry is built on the CPU for the whole grid and partitioned across
  ranks. Only the connected world ocean is kept, so the Black and Caspian Seas are land. `connect_basins!` keeps
  the Mediterranean: if the Gulf of Cádiz and the Alboran Sea fall in different basins, it deepens to the 280 m
  Camarinal Sill the cells along the least-deepening, face-connected path between them. In practice the strait is
  closed at 1°, 1/3° and 1/6° (2–4 cells deepened) and already open at 1/2°, 1/4° and 1/12°.
  Dead ends are filled: each column is raised to the depth of its second-deepest neighbor until every cell has at
  least two open sides, since seafloor-absorbed shortwave heated one-sided cells (pits, one-cell inlets) to 35–37 ᵒC.
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
- **Inputs**: each run job first stages what it reads onto `/work/nvme` (`slurm/stage_inputs.sh`): the GLORYS initial
  condition, 2015 JRA55 runoff, the bathymetry cache, and the segment's ERA5 hours, with symlinks to `/work/hdd` for
  everything else (about 20 GB). Cold reads from `/work/hdd` made the model build take up to 27 minutes.
- **Output** (`$RUN_DIRECTORY/<name>/`, on `/work/nvme`): daily surface T, S, e, u, v, w; η; sea ice h, ℵ, u, v. Checkpoints every
  5 simulated days, so 1-hour interactive jobs of 5 days each can be chained. Each run picks up from the latest checkpoint, so you extend a run by resubmitting with a
  larger `stop_days`.

## Pinned versions

Julia ≥ 1.12.3 is required: earlier versions assume 4 KiB pages on aarch64, and JLD2's memory-mapped writes then
fail on DeltaAI's 64 KiB pages (`SystemError("msync", 22)`, JLD2 #702).

`Project.toml` pins NumericalEarth at `3a1ee30` (branch `glw/global-kilo-mischief`), which is `main` plus eight open PRs:

- #742: fixes GLORYS inpainting. Deep levels and land were left at zero, and atoll columns got surface water copied down to the seafloor.
- #743: lets ERA5 supply mean sea level pressure instead of surface pressure.
- #769: builds JRA55 river routing on distributed grids (it crashed rebuilding the grid on the CPU).
- #770: adds `frazil_formation_depth`. We set 50 m: deep cells cooled below freezing by advection undershoots at steep
  topography (900 m near the Faroes) made up to 5.5 m/day of surface ice over 9 ᵒC water.
- #773: reads each ERA5 window's files concurrently before loading them, and the next window's in the background,
  cutting a 48-hour reload from about 70 s in a run.
- #774: searches river-mouth target cells within a latitude band, cutting `JRA55PrescribedLand` at 1/6° from 364 s to 31 s.
- #776: reads each ERA5 variable once when splitting a CDS response into hourly files (0.96 s → 0.03 s per file).
- #779: gives split ERA5 files their final name only once complete, so a killed download cannot leave a corrupt file.

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

sbatch slurm/download_data.sbatch 180                        # bathymetry + forcing for 180 days
sbatch slurm/smoke_test.sbatch
sbatch slurm/sixth_degree.sbatch 10 20                       # stop_days max_Δt_minutes; rerun with more days to extend
sbatch slurm/twelfth_degree.sbatch 10 10
```

## Known limitations

- JRA55-do runoff ends on 2019-12-31, which caps hindcasts started in 2015 at five years.
- No tides yet. They are planned once the runs are stable, as an equilibrium tidal potential through the
  pluggable barotropic potential of NumericalEarth #736.
