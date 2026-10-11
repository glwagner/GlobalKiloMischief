# Progress

## Status

| | 1/6ᵒ, 1 GPU | 1/12ᵒ, 4 GPUs |
|---|---|---|
| Furthest stable run | 30 days at Δt = 20 min | not yet run on GPUs |
| SYPD | 2.05 between ERA5 reloads, about 1.5 including them | queued (1-hour interactive test) |
| Largest Δt tried | 20 min (cap, not a stability limit) | 10 min planned |

**Open issues**
- Bottom cells in single-cell bathymetry pits heat to 37 ᵒC (Amazon mouth).
- Sea ice reaches 11–13 m in a few enclosed coastal cells.
- ERA5 reload costs about 27% of throughput.

## Log

**2026-10-11**
- Day 30 reached at 1/6ᵒ, stable.
- Warm maximum (37.1 ᵒC): bottom cells of single-cell pits, whose four horizontal neighbors are solid, at the Amazon
  mouth. Next: fill pits in `global_bottom_height`.
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
