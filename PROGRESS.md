# Progress

## Status

| | 1/6ᵒ, 1 GPU | 1/12ᵒ, 4 GPUs |
|---|---|---|
| Furthest stable run | day 140 (Δt ≤ 20 min to day 35, ≤ 30 to 45, ≤ 45 to 55, ≤ 45 at CFL 0.7 after) | 1 day at Δt = 10 min |
| SYPD | Δt ≤ 20 min: 1.95. ≤ 30: 2.5. ≤ 45: 3.1. ≤ 45 at CFL 0.7: **3.4** (Δt 25–45 min, 1.78 s/step) | 0.74 (2.2 s/step), first 100 steps |
| Largest Δt tried | 45 min cap: stable 10 days at CFL 0.5 (Δt ≈ 35 min); 4 days at CFL 0.7 (Δt ≈ 44 min) | 10 min, from rest (cap) |

**Open issues**
- Sea ice reaches 11–13 m in a few enclosed coastal cells.
- Per-step overheads: ERA5 reload 5–11 s per 47 h (was 70 s), daily output 1.8 s, checkpoint 24 s, and a one-time
  ~160 s stall at the first daily output after a pickup (output writers write the grid on first use).
- A 1-hour job spends about 17 min before the first step: model build 10.4 min (ERA5 first window 2–4.5 min, compilation ~5), grid 1.2, first step 5.6.

## Log

**2026-10-15**
- Days 110–140 stable (three segments): Δt 23–45 min, T in (−2.03, 35.7) ᵒC, sea ice back to ~16 m.
- Segment ends now fall on the 5-day checkpoints, so no simulated days are redone.
- A download failed opening a corrupt file inside the CDS zip (HDF error): transient, retried.

**2026-10-14**
- Days 100–114 stable at Δt ≈ 45 min. Tmax 35.8 ᵒC, sea ice up to 25.5 m. Day-110 checkpoint shows the warm
  cells are shallow (10–20 m) dead-end bottom cells, walled on three sides, in Sumatra (35.6 ᵒC), the Gulf of
  Martaban (35.1 ᵒC) and Flores: little exchange plus bottom-absorbed shortwave, like the pits. Possible fix: raise
  each column to its second-deepest neighbor, so every cell has two open sides. Not applied (changes the grid).
  Thickest ice: one 20 m coastal cell at the Filchner front (−40ᵒ, −78ᵒ).
- Download job segfaulted in `nc_close` after ~2000 files; its `.partial` file was discarded, no corrupt files.
- Days 85–104 stable, mostly at Δt = 45 min (1.7–1.9 s/step); T in (−2.34, 35.1) ᵒC; sea ice up to 17.3 m.
  Stopped at day 104 on another file a killed download left (u10 2015-04-15T03: full size, data unreadable).
  Fix: [NumericalEarth #779](https://github.com/NumericalEarth/NumericalEarth.jl/pull/779) writes `.partial` then renames; pinned `3a1ee30`.
- CFL 0.7, days 70–85: stable, T in (−2.47, 34.9) ᵒC; sea ice spiked to 18.6 m for a day, back to 13.5 m. Δt 27–45 min.
- Forcing downloads and run segments now alternate in the interactive slot, chained by `logs/download_then_run.sh`.

**2026-10-13**
- CFL 0.7, days 55–74: stable, T in (−2.28, 34.8) ᵒC, sea ice ≤ 14.1 m. Δt averages ~37 min, dipping to 25 min when
  max |u| reaches 4–5 m/s. ~3.4 SYPD. Stopped at day 74 on a truncated ERA5 file (`strd`, 2015-03-17T22) that a killed
  download job had left. Removed it and an empty one; checkpoint at day 70.
- Fix: split files now take their final name only when complete (NumericalEarth branch `glw/era5-atomic-split`, testing).
- CFL 0.7, max Δt 45 min, days 55–59: stable, Δt 41–45 min, ~1.74 s/step → **~4.1 SYPD** (short test, few samples).
- Startup is growing with `/work/hdd` contention: model build 17 min (ERA5 first window 3 + 3 min), first step 8.8 min.
- ERA5 split is 100× too slow: CDS chunks span 100 hours and are decompressed once per hourly file (0.96 s vs
  0.01 s per hour read in full). Fixed: 0.03 s per file, tests pass
  ([NumericalEarth #776](https://github.com/NumericalEarth/NumericalEarth.jl/pull/776), pinned `82858a6`). Download to day 180 rerunning.

**2026-10-12**
- Fastest currents (day 55): narrow shallow straits one or two cells wide: Bering 2.7 m/s, Karimata (Indonesia) 2.5,
  Belle Isle 2.2, Skagerrak 1.6. Gulf Stream and Florida Current 1.7–1.8 m/s are realistic. These set the CFL limit.
  Next: test `cfl = 0.7` (now a script argument) with max Δt 45 min.
- ERA5 was cached only to day 61. The login-node download hit 15.5 GB (CDS sends ~24 days of all variables as one
  7 GB file, split in memory) and was killed. Compute nodes reach CDS, but the split writes ~1 hourly file/s, so 24 days
  take ~90 min: now a 4-hour `ghx4` job (`slurm/download_data.sbatch`). Meanwhile, CFL 0.7 test over days 55–59.
- Δt ≤ 45 min, days 45–55: stable, T in (−2.14, 34.6) ᵒC, sea ice up to 14.6 m. The CFL limit keeps Δt at
  29–38 min; 1.7–1.8 s/step (more barotropic substeps) → 3.1 SYPD. Continuing this as the main run.
- ERA5 first-window load at startup varies 3–8 min (cold `/work/hdd`).
- Δt ≤ 30 min, days 35–45: stable, T in (−2.14, 34.6) ᵒC. 1.64 s/step at 30 min = 3.0 SYPD. The CFL wizard
  drops Δt to 19–23 min when coastal currents reach 4.2 m/s, so the average is about 2.5 SYPD. Next: Δt ≤ 45 min.
- 1/6ᵒ days 25–35 in 43 min, stable: T in (−2.34, 34.5) ᵒC, sea ice ≤ 13.5 m, falling to 10.6 m. With #773 and #774:
  JRA55 land 364 → 29 s, ERA5 reloads 33 → 5–11 s after the first. **1.95 SYPD sustained.**
- Branched a max Δt = 30 min run from the day-35 checkpoint (`sixth_degree_dt30min`), 10 days.
- 1/6ᵒ days 15–25, stable: T in (−2.23, 33.8) ᵒC, sea ice ≤ 12.7 m. Slow-step log attributes the overheads (status above).
- JRA55 land build: 364 s → 31 s. Brute-force nearest-wet-cell search for 16k river mouths; now searched within a
  latitude band, results identical ([NumericalEarth #774](https://github.com/NumericalEarth/NumericalEarth.jl/pull/774)).
  The in-run "JRA55 stall" was the output writers' first write, not JRA55 (its reloads take 0.1 s).
- #773 now also reads the next ERA5 window's files in the background after each reload. Both pinned (`77de1e0`).
- 1/6ᵒ days 5–15 in 48 min, stable: T in (−2.03, 33.3) ᵒC, sea ice up to 13 m (fold/coastal cells again).
  1.62 SYPD over the segment, up from about 1.5. Step times still spike irregularly (one 50-step window cost +125 s, others +25 s).
  Next segment reports each slow step and times each stage of the model build.
- 1/12ᵒ runs on 4 GPUs: 1 day, stable, 2.2 s/step at Δt = 10 min, **0.74 SYPD** (72% weak-scaling efficiency vs 1/6ᵒ).
  T in (−2.0, 34.8) ᵒC, max|u| 4.2 m/s, 78 GiB per GPU. First step compiles for 19 min.
- It started at Δt = 10 min from rest, by accident (it treated the iteration-0 checkpoint as a pickup), and held.
  Fixed: only checkpoints past iteration 0 skip the ramp.
- 1/6ᵒ restart with filled pits: day 8 in 1 hour, stable. The warm maximum falls from 35.6 ᵒC to 33.2 ᵒC (it rose
  to 36.8 ᵒC with pits). Minimum −2.03 ᵒC.
- The ERA5 prefetch (#773) never ran: the manifest kept the previous commit's tree hash, so Julia loaded the old
  NumericalEarth. Reload steps still cost about 70 s. Reinstalled; next segment measures it.

**2026-10-11**
- 1/12ᵒ test failed at the first time step after 29 min of setup: the 4 ranks created different numbers of
  distributed fields. Cause: in the `progress` callback, `maximum(h * ℵ)` builds one extra field on the fold rank of a
  distributed tripolar grid (an Oceananigans bug that affects any reduction of a 2D field). Workaround: take the
  maximum of the interior array and combine it with `all_reduce`; checked on a 4-rank toy grid. Resubmitted the 1-hour
  4-GPU test.
- ERA5 reload: almost all of the time is cold reads from `/work/hdd`. Reading a window's files concurrently cuts a
  reload from about 50 s to 15 s ([NumericalEarth #773](https://github.com/NumericalEarth/NumericalEarth.jl/pull/773),
  now pinned). Julia threads now follow the Slurm CPU count.
- Pits filled: every column deeper than all four neighbors is raised to its deepest neighbor (37k columns at 1/6ᵒ,
  63k at 1/12ᵒ). Shortwave absorbed at the seafloor had heated pits at the Amazon mouth to 37 ᵒC. The 1/6ᵒ run
  restarts from day 0 after the 1/12ᵒ test.
- Day 30 reached at 1/6ᵒ, stable.
- Sea-ice maximum (11.3 m): an Antarctic coastal cell 10 m deep, mostly enclosed by land.
- 1/12ᵒ: submitted a 1-hour 4-GPU interactive test. Job 3354360 on `ghx4` stays queued (estimated start 2026-10-12).
- Picked-up runs now start at `max_Δt`. Checkpoints do not store Δt, so each run was ramping up from 1 min again.
- Checkpoints every 5 days, so 10-day segments fit in 1-hour interactive jobs, which start within minutes.

**2026-10-10**
- Cold cells (−12.7 ᵒC) in the Faroe Bank Channel and Denmark Strait overflows. Cause: next to walls, tracer WENO
  falls back to second-order centered reconstruction, which leaves bottom corner cells undamped. Fix:
  `minimum_buffer_upwind_order = 1` (upwind fallback). The minimum now stays near −2 ᵒC. Upstream:
  [NumericalEarth #772](https://github.com/NumericalEarth/NumericalEarth.jl/pull/772).
- The same cause made the spurious Faroes ice. Frazil formation is now limited to the upper 50 m (NumericalEarth #770).
- Sea ice uses incremental remapping instead of WENO, which removed a 69 m pile-up at Adélie Land.
- 1/6ᵒ throughput: 1.56–1.60 s per step at Δt = 20 min.
