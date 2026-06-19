# Observable-driven flux map + indexed explorer lookups

**Date:** 2026-06-19
**Status:** Approved (design)

## Problem

The interactive explorer (`explore()`) feels unresponsive even though all ODE
solves are precomputed and cached as CSVs. The grids are tiny (the 3-D
`atp_demand_grid.csv` is ~1,584 rows), so lookups are microseconds — this is a
**rendering** problem, not an optimization/solver problem.

The root cause is in the explorer's `onany(atp_level, selected)` callback
(`src/explorer.jl`). On *every* ATP-slider tick and *every* cell click it runs:

```julia
empty!(ax_net)
draw_ppp_flux_map!(ax_net, fluxes; title = ttl)
```

`empty!` tears down the entire network scene and `draw_ppp_flux_map!` rebuilds
it from scratch. That rebuild is not a few plots — it is ~80–100 individual
primitives created in loops (`src/flux_map.jl`):

- per reaction leg: a `lines!` (shaft) + a `scatter!` (arrowhead)
- per reaction: a `poly!` + `text!` (enzyme square)
- one `scatter!` **per node** (~17, not batched) + one `text!` per node label
- background group `poly!`s

Dragging the slider fires this callback continuously, so the scene is fully
destroyed and rebuilt many times per second. That teardown/rebuild is the jank.

Secondary (negligible but on the same hot path): `zmatrix` rescans all ~1,584
rows on every slider tick, and the callback's `filter(eachrow(df))` is a full
scan to fetch one row.

## Goal

Make the explorer responsive by building the network plot objects **once** and
driving them with `Observable`s (the idiomatic Makie pattern), so an interaction
just pushes new values into existing plots — no teardown/rebuild. Additionally
replace the per-interaction dataframe scans with pre-built indices.

The static `phasemap()` path and all existing tests must keep working unchanged.

## In-repo precedent

`src/metab_panel.jl` is **already** Observable-driven: `draw_metab_panel!(ax,
lit, metab_obs)` draws the static literature boxes once and positions the dynamic
model diamonds via `lift(metab_obs)`. The explorer's metab panel therefore
already updates the right way; only the flux map fights it. This refactor makes
the flux map match that pattern.

## Design

### Component 1 — `flux_map.jl`: Observable-driven core

Split `draw_ppp_flux_map!` into static-once structure and dynamic-via-`lift`
geometry.

- **New core signature:**
  `draw_ppp_flux_map!(ax, fluxes_obs::Observable; title = "", maxflux = nothing)`.
- **Drawn once** (independent of flux value): background group polys, enzyme
  squares + labels, shared-pool halos, node dots, node labels.
- **Dynamic:** a single top-level
  `legdata = lift(fluxes_obs) do f ... end`
  computes `gmax`/`pmax` once per update and returns a **fixed-length** `Vector`
  of per-leg draw specs `(shaft_pts::Vector{Point2f}, lw, tip, color,
  head_pt::Point2f, head_rot)`, ordered by a flat enumeration of
  `(reaction, leg)`. Because `_PPP_REACTIONS` is a compile-time constant, the leg
  count is fixed, so exactly one `lines!` (shaft) + one `scatter!` (arrowhead)
  are created **per leg at build time**, each attribute `lift`ed off
  `legdata[j]`. An update recomputes the lifts; no plot is ever added or removed.
- **Geometry helper:** extract the current inline straight/curved + gap-trim +
  sign-flip logic (the body of the `for r in _PPP_REACTIONS` leg loop and
  `_draw_curved_flux!`/`_draw_flux_leg!` math) into a pure `_leg_drawspec(...)`
  returning one leg's spec, so the static and dynamic paths share identical math.
- **Backward-compat wrapper:**
  `draw_ppp_flux_map!(ax, fluxes::NamedTuple; title, maxflux) =
   draw_ppp_flux_map!(ax, Observable(fluxes); title, maxflux)`.
  Under CairoMakie each `lift` evaluates once, so `src/phasemap.jl` and
  `test/test_flux_map.jl` are untouched and there is one source of truth.
- **Title:** accept `title` as `String` *or* `Observable{String}`. When it is an
  Observable, bind with `on(title; update = true) do t; ax.title = t end`;
  otherwise set `ax.title` once.

### Component 2 — `explorer.jl`: build-once + indexed lookups

- **Row index, built once:** `Dict{Tuple{Float64,Int,Int}, DataFrameRow}` keyed
  `(atpase_frac, i_nadph, i_r5p)`. Exact `==` on the CSV-parsed `Float64`s is
  already documented-safe in this file. Replaces the per-interaction
  `filter(eachrow(df))`.
- **Z matrices, precomputed once:** `Dict(atp => Z)` (~11 ATP levels × 12×12
  grid). The `zmatrix` `lift` becomes a dict lookup rather than a full rescan.
- **Network built once:** create `fluxes_obs = Observable(...)` and
  `title_obs = Observable("")`, then call
  `draw_ppp_flux_map!(ax_net, fluxes_obs; title = title_obs)` **before** the
  callback.
- **Slimmed `onany(atp_level, selected)` callback:** index-lookup the row, then
  `fluxes_obs[] = row_to_fluxes(row)`, `title_obs[] = ttl`, set `footer[]`, set
  `metab_obs[]`. No `empty!`, no redraw.

## Data flow (after refactor)

```
slider / click
  → onany updates selected/atp_level Observables
    → row = index[(atp, i_nadph, i_r5p)]        # O(1) dict lookup
    → fluxes_obs[] = row_to_fluxes(row)         # push only
       → legdata lift recomputes gmax/pmax + per-leg specs
          → each leg's shaft/head lines!/scatter! lift updates in place
    → title_obs[], footer[], metab_obs[] pushed
  → heatmap Z_obs lift = zdict[atp]             # O(1) dict lookup
```

## Testing

- **Unchanged paths must still pass:** `test/test_flux_map.jl` (NamedTuple
  signature → wrapper) and the `phasemap` tests. This proves the wrapper and
  CairoMakie single-evaluation behavior.
- **Regression guard in `test/test_explorer.jl`:** after `build_explorer`,
  capture `length(ax_net.scene.plots)`; push a new `atp_level` and a new
  `selected`; assert the plot count is **unchanged** (proves built-once, not
  rebuilt) and that no error is thrown.
- **Index unit tests:** the `(atp, i_nadph, i_r5p)` lookup returns the same row
  the old `filter` did; the precomputed Z dict matches the old `zmatrix` output
  for every ATP level.

## Risks

The only real risk is Makie attribute-binding correctness: scalar
`linewidth`/`color`/`rotation` and a `Vector{Point2f}` shaft each driven by an
Observable. `metab_panel.jl` already proves `lift`→`scatter!` works in this
codebase; `lines!` with an Observable point vector plus scalar Observable
attributes is standard Makie. The plot-count regression test catches any
accidental rebuild.

## Out of scope

- No changes to the kinetic model, grid contents, or the extension.
- No slider throttling — Observable updates are cheap enough that it is
  unnecessary.
