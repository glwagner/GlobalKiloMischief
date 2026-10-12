# Progress

## Status

| | 1/6ᵒ, 1 GPU | 1/12ᵒ, 4 GPUs |
|---|---|---|
| Furthest stable run | day 15 of the restart with filled pits (30 days before) | 1 day at Δt = 10 min |
| SYPD | 2.1 between reloads, 1.62 over days 5–15 | 0.74 (2.2 s/step), first 100 steps |
| Largest Δt tried | 20 min (cap, not a stability limit) | 10 min, from rest (cap) |

**Open issues**
- Sea ice reaches 11–13 m in a few enclosed coastal cells.
- ERA5 reload costs about 25% of throughput (fix now installed, not yet measured in a run).
- A 1-hour job spends about 21 min before the first step: building the model takes 14 min (61% of it not compiling), the first step 5.6 min.

## Log

**2026-10-12**
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
