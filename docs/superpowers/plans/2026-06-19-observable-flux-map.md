# Observable-Driven Flux Map + Indexed Explorer Lookups — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the interactive explorer responsive by building the network flux-map plot objects once and driving them with Observables (instead of `empty!`/full-rebuild on every interaction), and replace the per-interaction dataframe scans with pre-built indices.

**Architecture:** Split `draw_ppp_flux_map!` into static-once structure plus dynamic per-leg geometry driven by a single `lift` off a `fluxes` Observable (mirroring the existing `draw_metab_panel!` pattern). Keep the `NamedTuple` signature as a thin wrapper so the static `phasemap()` path and existing tests are untouched. In the explorer, build the network once against a `fluxes_obs`/`title_obs`, slim the `onany` callback to pushes only, and serve row + Z lookups from `Dict` indices built once.

**Tech Stack:** Julia 1.10, Makie (GLMakie interactive / CairoMakie offscreen), DataFrames, CSV.

## Global Constraints

- Julia 1.10 is the minimum/tested version.
- Core `[deps]` stay exactly `CSV`, `DataFrames`, `Makie`, `GLMakie`, `CairoMakie` — **add no dependencies**.
- The viewer stays PPP-free: no `import`/`using` of `PentosePhosphatePathway` or any weakdep in core source.
- Tests render **offscreen with CairoMakie** (no GL window); they must run headless.
- No changes to the kinetic model, the package extension, or the cached grid CSVs.
- Follow each source file's existing comment style (`# file.jl — …` headers); do not introduce new banner styles.

---

### Task 1: Observable-driven `draw_ppp_flux_map!` core + compat wrapper

**Files:**
- Modify: `src/flux_map.jl` (replace `_draw_flux_leg!`/`_draw_curved_flux!` usage and the `draw_ppp_flux_map!` body; lines ~98–230)
- Test: `test/test_flux_map.jl`

**Interfaces:**
- Consumes: `_PPP_NODES`, `_PPP_REACTIONS`, `_PPP_GROUPS`, `_PPP_SHARED`, `_PPP_LABELS`, `_PPP_LABEL_OFF`, `_LEG_GAP`, `_draw_enzyme_square!` (all existing in `flux_map.jl`); `Observable`, `lift`, `on` (in scope via `using Makie`).
- Produces:
  - `_curved_pts(A, B, M, gA, gB) -> Vector{Tuple{Float64,Float64}}` — trimmed Bézier polyline points.
  - `_leg_drawspec(from, to, ctrl, v, lw, tip) -> NamedTuple` with fields `shaft::Vector{Point2f}`, `col::Symbol`, `lw::Float64`, `tip::Float64`, `head::Point2f`, `headrot::Float64`.
  - `const _PPP_LEGS::Vector{NamedTuple}` — fixed flat list of legs, fields `key`, `region`, `from`, `to`, `ctrl`.
  - `draw_ppp_flux_map!(ax, fluxes_obs::Observable; title = "", maxflux = nothing) -> ax` (core).
  - `draw_ppp_flux_map!(ax, fluxes::NamedTuple; title = "", maxflux = nothing) -> ax` (wrapper).

- [ ] **Step 1: Write the failing test**

Add this block at the END of the existing `@testset "draw_ppp_flux_map!" begin … end` in `test/test_flux_map.jl`, just before its closing `end` (so it reuses `base`/`rev` already defined in that testset):

```julia
    # Observable path: plots are built once and updated in place (no rebuild on push)
    obs = Observable(base)
    fig3 = Figure(size = (500, 460)); ax3 = Axis(fig3[1, 1])
    @test_nowarn G.draw_ppp_flux_map!(ax3, obs; title = "obs")
    plots_before = copy(ax3.scene.plots)
    @test_nowarn (obs[] = rev)
    @test length(ax3.scene.plots) == length(plots_before)
    @test all(a === b for (a, b) in zip(ax3.scene.plots, plots_before))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'using Pkg; Pkg.activate("."); include("test/test_flux_map.jl")'`
Expected: FAIL — the current untyped `draw_ppp_flux_map!(ax, fluxes; …)` runs `getfield` on the `Observable`, erroring (no `Observable` method yet).

- [ ] **Step 3: Write minimal implementation**

In `src/flux_map.jl`, **delete** the two helpers `_draw_flux_leg!` (currently ~lines 106–113) and `_draw_curved_flux!` (currently ~lines 121–133), and **replace** the entire `draw_ppp_flux_map!` function (currently ~lines 146–230) with the following. Keep everything above `_draw_flux_leg!` (the node/reaction constants and `_LEG_GAP`) and the `_draw_enzyme_square!` helper unchanged.

```julia
# Trimmed cubic Bézier polyline from node centre A (tail) to B (head) bowing through M (the
# enzyme square centre), with `gA`/`gB` data-unit gaps trimmed off each end so the shaft starts/
# stops shy of the node dots. Returns the point list (tuples). The handle length s = 4/3 makes
# the curve pass through M at t = 0.5 so paired bi-bi legs meet behind the square.
function _curved_pts(A, B, M, gA, gB)
    s = 4 / 3
    P1 = (A[1] + s * (M[1] - A[1]), A[2] + s * (M[2] - A[2]))
    P2 = (B[1] + s * (M[1] - B[1]), B[2] + s * (M[2] - B[2]))
    bez(t) = (u = 1 - t; (u^3*A[1] + 3u^2*t*P1[1] + 3u*t^2*P2[1] + t^3*B[1],
                          u^3*A[2] + 3u^2*t*P1[2] + 3u*t^2*P2[2] + t^3*B[2]))
    pts = [bez(t) for t in range(0, 1; length = 64)]
    i0 = findfirst(p -> hypot(p[1] - A[1], p[2] - A[2]) >= gA, pts)
    i1 = findlast(p  -> hypot(p[1] - B[1], p[2] - B[2]) >= gB, pts)
    i0 === nothing && (i0 = 1)
    (i1 === nothing || i1 <= i0) && (i1 = length(pts))
    return pts[i0:i1]
end

# Pure geometry for ONE flux leg. Returns the shaft polyline (as Point2f), colour (reverse =
# steel blue, forward = tomato), width/tip pixel sizes, and the arrowhead position + rotation.
# A negative flux flips tail<->head so the arrowhead points the reverse way. LOWER is a text-only
# label (no dot), so its end takes no gap.
function _leg_drawspec(from, to, ctrl, v, lw, tip)
    col = v < 0 ? :steelblue : :tomato
    nf, nt = from, to
    A = _PPP_NODES[from]; B = _PPP_NODES[to]
    if v < 0
        (A, B) = (B, A); (nf, nt) = (nt, nf)
    end
    gA = nf === :LOWER ? 0.0 : _LEG_GAP
    gB = nt === :LOWER ? 0.0 : _LEG_GAP
    if ctrl === nothing
        L = hypot(B[1] - A[1], B[2] - A[2])
        if L > gA + gB
            ux, uy = (B[1] - A[1]) / L, (B[2] - A[2]) / L
            A = (A[1] + ux * gA, A[2] + uy * gA)
            B = (B[1] - ux * gB, B[2] - uy * gB)
        end
        pts = [A, B]
    else
        pts = _curved_pts(A, B, ctrl, gA, gB)
    end
    he = pts[end]; ref = pts[max(1, lastindex(pts) - 1)]
    headrot = atan(he[2] - ref[2], he[1] - ref[1]) - pi/2   # utriangle points +y at rotation 0
    shaft = Point2f[Point2f(p[1], p[2]) for p in pts]
    return (shaft = shaft, col = col, lw = lw, tip = tip,
            head = Point2f(he[1], he[2]), headrot = headrot)
end

# Fixed flat list of legs (built once at load): the leg COUNT never changes, so the explorer can
# build exactly one shaft + one arrowhead plot per entry and only update their data on reselect.
const _PPP_LEGS = [(key = r.key, region = r.region, from = leg[1], to = leg[2], ctrl = leg[3])
                   for r in _PPP_REACTIONS for leg in r.legs]

"""
    draw_ppp_flux_map!(ax, fluxes_obs::Observable; title = "", maxflux = nothing)
    draw_ppp_flux_map!(ax, fluxes::NamedTuple; title = "", maxflux = nothing)

Draw the glycolysis+PPP wiring. The static structure (region tints, enzyme squares, shared-pool
halos, node dots, labels) is drawn once. The reaction legs are driven by `fluxes_obs`: a single
`lift` recomputes per-leg geometry (colour = sign, width ∝ |flux|, arrowhead direction = sign) and
each leg's shaft (`lines!`) + arrowhead (`scatter!`) updates IN PLACE — no plot is added or removed
on reselect. Glycolysis and PPP self-normalize on separate width scales (PPP max width = ½
glycolysis) unless `maxflux` overrides the glycolysis reference. `title` may be a String or an
`Observable{String}`. The `NamedTuple` method wraps the value in a constant `Observable`, so the
static CairoMakie path evaluates each `lift` exactly once.
"""
function draw_ppp_flux_map!(ax, fluxes_obs::Observable; title = "", maxflux = nothing)
    hidedecorations!(ax); hidespines!(ax)
    if title isa Observable
        on(title; update = true) do t
            ax.title = t
        end
    else
        ax.title = title
    end
    xlims!(ax, 0.0, 7.1); ylims!(ax, 0.1, 5.7)

    # static background region tints
    for (xlo, xhi, ylo, yhi, col) in _PPP_GROUPS
        poly!(ax, Point2f[(xlo, ylo), (xhi, ylo), (xhi, yhi), (xlo, yhi)]; color = col)
    end

    # dynamic per-leg draw specs: ONE lift; gmax/pmax computed once per update. Two independent
    # width scales so the PPP branch stays legible (glycolysis max 14px; PPP max 7px = ½).
    legdata = lift(fluxes_obs) do f
        _gm(pred) = begin
            m = maximum((abs(getfield(f, r.key)) for r in _PPP_REACTIONS if pred(r.region));
                        init = 0.0) do v
                isnan(v) ? 0.0 : v
            end
            (m == 0 || isnan(m)) ? 1.0 : m
        end
        gmax = maxflux === nothing ? _gm(==(:glyc)) : maxflux
        gmax = (gmax == 0 || isnan(gmax)) ? 1.0 : gmax
        pmax = _gm(!=(:glyc))
        map(_PPP_LEGS) do L
            v = getfield(f, L.key)
            av = isnan(v) ? 0.0 : abs(v)
            lw = L.region === :glyc ? 1.0 + 13.0 * av / gmax : 1.0 + 6.0 * av / pmax
            tip = 8.0 + 2.2 * lw
            _leg_drawspec(L.from, L.to, L.ctrl, v, lw, tip)
        end
    end

    # one shaft (lines!) + one arrowhead (scatter!) per leg, built ONCE, lifted by index. Drawn
    # AFTER the region tints and BEFORE the squares so legs stay behind the enzyme squares.
    for j in eachindex(_PPP_LEGS)
        col = lift(d -> d[j].col, legdata)
        lines!(ax, lift(d -> d[j].shaft, legdata); color = col,
               linewidth = lift(d -> d[j].lw, legdata))
        scatter!(ax, lift(d -> d[j].head, legdata); marker = :utriangle,
                 markersize = lift(d -> d[j].tip, legdata),
                 rotation = lift(d -> d[j].headrot, legdata), color = col)
    end

    # static enzyme squares on top (+ any `extra` label, e.g. a reaction's separated arm)
    for r in _PPP_REACTIONS
        _draw_enzyme_square!(ax, r.square, r.label, r.region)
        ex = get(r, :extra, nothing)
        ex === nothing || _draw_enzyme_square!(ax, ex[1], ex[2], r.region)
    end
    # static shared-pool halos, then nodes (LOWER is text-only), then labels
    for name in _PPP_SHARED
        (x, y) = _PPP_NODES[name]
        scatter!(ax, [x], [y]; markersize = 26, color = (:gold, 0.25))
    end
    for (name, (x, y)) in _PPP_NODES
        name === :LOWER && continue
        scatter!(ax, [x], [y]; markersize = 12, color = :white,
                 strokecolor = :black, strokewidth = 1.8)
    end
    for (name, (x, y)) in _PPP_NODES
        off = get(_PPP_LABEL_OFF, name, (0.0, 0.42, (:center, :bottom)))
        text!(ax, x + off[1], y + off[2]; text = get(_PPP_LABELS, name, string(name)),
              align = off[3], fontsize = 15,
              font = (name in _PPP_SHARED ? :bold : :regular))
    end
    return ax
end

draw_ppp_flux_map!(ax, fluxes::NamedTuple; title = "", maxflux = nothing) =
    draw_ppp_flux_map!(ax, Observable(fluxes); title = title, maxflux = maxflux)
```

- [ ] **Step 4: Run the flux-map test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.activate("."); include("test/test_flux_map.jl")'`
Expected: PASS (both the existing NamedTuple smoke tests and the new Observable in-place test).

- [ ] **Step 5: Run the full suite to confirm the static `phasemap` path is unaffected**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS — `phasemap`/`render_phase_map!` tests still green (they hit the `NamedTuple` wrapper; CairoMakie evaluates each `lift` once).

- [ ] **Step 6: Commit**

```bash
git add src/flux_map.jl test/test_flux_map.jl
git commit -m "$(printf 'refactor(flux_map): Observable-driven core, build-once legs\n\nSplit draw_ppp_flux_map! into static structure (drawn once) and dynamic\nper-leg geometry driven by a single lift off a fluxes Observable, so the\nlegs update in place instead of being rebuilt. Keep the NamedTuple\nsignature as a constant-Observable wrapper so phasemap and existing tests\nare unchanged.\n\nCo-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>')"
```

---

### Task 2: Explorer lookup indices (`_build_row_index`, `_build_zdict`)

**Files:**
- Modify: `src/explorer.jl` (add two helper functions near the top, after the `logclick_to_cell` helper ~line 6)
- Test: `test/test_explorer.jl`

**Interfaces:**
- Consumes: a loaded grid `df` (columns `atpase_frac`, `i_nadph`, `i_r5p`, `retcode::Symbol`, `cycle_index`).
- Produces:
  - `_build_row_index(df) -> Dict{Tuple{Float64,Int,Int}, <:DataFrameRow}` keyed `(atpase_frac, i_nadph, i_r5p)`.
  - `_build_zdict(df, nx, ny) -> Dict{Float64, Matrix{Float64}}` — per ATP level, the `[i_r5p, i_nadph]` cycle-index matrix; masked/non-converged cells are `NaN`.

- [ ] **Step 1: Write the failing test**

Add this new testset at the END of `test/test_explorer.jl` (after the existing `@testset "build_explorer offscreen"` block):

```julia
@testset "explorer lookup indices" begin
    df = G.load_grid(G.default_atp_grid())
    nx = maximum(df.i_r5p); ny = maximum(df.i_nadph)

    # row index returns the same row the old per-interaction filter did
    idx = G._build_row_index(df)
    r = first(filter(x -> x.retcode == :Terminated, eachrow(df)))
    got = idx[(r.atpase_frac, r.i_nadph, r.i_r5p)]
    @test got.cycle_index == r.cycle_index
    @test got.r5p_phi == r.r5p_phi && got.nadph_phi == r.nadph_phi

    # zdict matches the old per-level zmatrix construction for EVERY ATP level
    zd = G._build_zdict(df, nx, ny)
    @test Set(keys(zd)) == Set(unique(df.atpase_frac))
    for atp in unique(df.atpase_frac)
        expected = fill(NaN, nx, ny)
        for rr in eachrow(df)
            (rr.atpase_frac == atp && rr.retcode == :Terminated) || continue
            expected[rr.i_r5p, rr.i_nadph] = rr.cycle_index
        end
        Z = zd[atp]
        @test size(Z) == (nx, ny)
        @test all((isnan.(Z) .& isnan.(expected)) .| (Z .== expected))
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'using Pkg; Pkg.activate("."); include("test/test_explorer.jl")'`
Expected: FAIL with `UndefVarError: _build_row_index not defined` (and `_build_zdict`).

- [ ] **Step 3: Write minimal implementation**

In `src/explorer.jl`, immediately AFTER the `logclick_to_cell(...)` one-liner (~line 6) and BEFORE the `build_explorer` docstring, insert:

```julia
# O(1) row lookup keyed (atpase_frac, i_nadph, i_r5p), built once. Exact == on the CSV-parsed
# Float64 atpase_frac is safe: the slider values come from the same parsed column.
_build_row_index(df) = Dict((r.atpase_frac, r.i_nadph, r.i_r5p) => r for r in eachrow(df))

# Per-ATP-level cycle-index matrices [i_r5p, i_nadph], built once. Masked/non-converged cells
# stay NaN (renderers mask on the :Terminated retcode). Replaces a full df rescan per slider tick.
function _build_zdict(df, nx, ny)
    zd = Dict(atp => fill(NaN, nx, ny) for atp in unique(df.atpase_frac))
    for r in eachrow(df)
        r.retcode == :Terminated || continue
        zd[r.atpase_frac][r.i_r5p, r.i_nadph] = r.cycle_index
    end
    return zd
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.activate("."); include("test/test_explorer.jl")'`
Expected: PASS — both the existing `build_explorer offscreen` testset and the new `explorer lookup indices` testset.

- [ ] **Step 5: Commit**

```bash
git add src/explorer.jl test/test_explorer.jl
git commit -m "$(printf 'feat(explorer): add row + Z lookup indices\n\nPre-build a (atpase_frac, i_nadph, i_r5p) row index and per-level\ncycle-index matrices so the explorer can serve reselect lookups in O(1)\ninstead of rescanning the full grid dataframe each interaction.\n\nCo-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>')"
```

---

### Task 3: Wire `build_explorer` to build the network once and push updates

**Files:**
- Modify: `src/explorer.jl` (`build_explorer`: replace the `zmatrix` closure + `Z_obs` lift, add the build-once network, slim the `onany` callback — currently ~lines 25–103)
- Test: `test/test_explorer.jl` (add a plot-identity regression guard to the existing `build_explorer offscreen` testset)

**Interfaces:**
- Consumes: `_build_row_index`, `_build_zdict` (Task 2); `draw_ppp_flux_map!(ax, fluxes_obs::Observable; title)` (Task 1); existing `row_to_fluxes`, `row_to_pools`, `_PPP_FLUX_KEYS`.
- Produces: unchanged public return `build_explorer(df) -> (fig, selected, atp_level)`; the network `Axis` is built once and updated in place.

- [ ] **Step 1: Write the failing test**

In `test/test_explorer.jl`, inside the existing `@testset "build_explorer offscreen" begin … end`, locate the two existing lines:

```julia
    @test_nowarn (selected[] = (conv.i_nadph, conv.i_r5p))
    @test_nowarn (atp_level[] = atp_levels[1])
```

Replace **just those two lines** with the following (captures the network axis's plots before mutating, then asserts the SAME plot objects persist after — the `empty!`/rebuild path would replace them, so this fails on current code):

```julia
    ax_net = only(filter(x -> x isa Axis && occursin("ATP", x.title[]), fig.content))
    net_plots_before = copy(ax_net.scene.plots)
    @test_nowarn (selected[] = (conv.i_nadph, conv.i_r5p))
    @test_nowarn (atp_level[] = atp_levels[1])
    @test length(ax_net.scene.plots) == length(net_plots_before)
    @test all(a === b for (a, b) in zip(ax_net.scene.plots, net_plots_before))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'using Pkg; Pkg.activate("."); include("test/test_explorer.jl")'`
Expected: FAIL on the `all(a === b …)` assertion — the current `onany` does `empty!(ax_net)` then redraws, so the post-update plots are NEW objects (identity differs).

- [ ] **Step 3: Write minimal implementation**

In `src/explorer.jl`, make three edits inside `build_explorer`:

**(3a)** DELETE the `zmatrix` closure (currently ~lines 25–32):

```julia
    # cycle-index Z matrix [r5p, nadph] for one ATP level; masked/non-converged cells stay NaN
    zmatrix(atp) = begin
        Z = fill(NaN, nx, ny)
        for r in eachrow(df)
            (r.atpase_frac == atp && r.retcode == :Terminated) || continue
            Z[r.i_r5p, r.i_nadph] = r.cycle_index
        end
        Z
    end
```

and in its place put the index construction:

```julia
    row_index = _build_row_index(df)
    zdict = _build_zdict(df, nx, ny)
```

**(3b)** Replace the `Z_obs` lift line (currently `Z_obs = lift(zmatrix, atp_level)`, ~line 48) with:

```julia
    Z_obs = lift(atp -> zdict[atp], atp_level)
```

**(3c)** Replace the build-once network setup and the whole `onany(...)` callback. The `selected` Observable and `marker_pos`/`scatter!` block (current ~lines 74–79) stay as-is. Immediately AFTER the `scatter!(ax_heat, marker_pos; …)` call and BEFORE the `onany` block, insert the build-once network; then replace the existing `onany(atp_level, selected) do atp, (in_, ir) … end` block (current ~lines 82–103) entirely. The combined replacement is:

```julia
    # build the network ONCE against Observables; reselect just pushes new values (no rebuild)
    fluxes_obs = Observable(row_to_fluxes(row_index[(atp_level[], selected[][1], selected[][2])]))
    title_obs = Observable("")
    draw_ppp_flux_map!(ax_net, fluxes_obs; title = title_obs)

    # redraw-free update of network + footer + metab on EITHER an ATP-level change or a cell click
    onany(atp_level, selected) do atp, (in_, ir)
        row = row_index[(atp, in_, ir)]
        atp_pct = round(atp * 100, sigdigits = 2)
        title_obs[] =
            "ATP $(atp_pct)%  |  NADPH demand $(round(row.nadph_phi * 100, sigdigits = 2))% / " *
            "R5P demand $(round(row.r5p_phi * 100, sigdigits = 2))% of supply"
        fluxes_obs[] = row_to_fluxes(row)
        mx = maximum((abs(getfield(fluxes_obs[], k)) for k in _PPP_FLUX_KEYS); init = 0.0) do v
            isnan(v) ? 0.0 : v
        end
        mx = (mx == 0 || isnan(mx)) ? 1.0 : mx
        footer[] = row.retcode == :Terminated ?
            "ATP $(atp_pct)%  |  cell (nadph=$in_, r5p=$ir)  |  mode=$(row.mode)  |  " *
            "cycle index=$(round(row.cycle_index, sigdigits = 3))  |  " *
            "max |flux| = $(round(mx, sigdigits = 3)) µM/min" :
            "ATP $(atp_pct)%  |  cell (nadph=$in_, r5p=$ir)  |  NON-CONVERGED ($(row.retcode))"
        metab_obs[] = row_to_pools(row)
    end
```

Note: the network build must come after `ax_net` (current ~line 58), `footer`/`metab_obs` (current ~lines 59, 68), and `selected` (current ~line 74) are defined — placing it right before `onany` satisfies all of these. The trailing `notify(selected)` (current ~line 112) stays and populates the title/footer/metab/fluxes for the initial cell.

- [ ] **Step 4: Run the explorer test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.activate("."); include("test/test_explorer.jl")'`
Expected: PASS — the plot-identity guard now holds (same objects, updated in place), and the index testset from Task 2 still passes.

- [ ] **Step 5: Run the full suite**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS across all test files.

- [ ] **Step 6: Regenerate the README screenshot to confirm `build_explorer` still renders end-to-end**

Run: `julia --project=. scripts/render_readme_screenshot.jl`
Expected: completes without error and rewrites the screenshot PNG. Review `git diff --stat` for the image; revert it if you don't intend to update the committed screenshot (`git checkout -- <png>`).

- [ ] **Step 7: Commit**

```bash
git add src/explorer.jl test/test_explorer.jl
git commit -m "$(printf 'perf(explorer): build network once, push updates in place\n\nReplace the per-interaction empty!/full-rebuild of the flux-map axis with\na build-once network driven by fluxes/title Observables, and serve row +\nZ lookups from the pre-built indices. The onany callback now only pushes\nnew values, so dragging the ATP slider no longer tears down and rebuilds\nthe scene each tick.\n\nCo-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>')"
```

---

## Self-Review

**Spec coverage:**
- Observable-driven core + static-once structure → Task 1 ✓
- Geometry helper shared by static/dynamic (`_leg_drawspec`/`_curved_pts`) → Task 1 ✓
- Fixed-count legs built once → Task 1 (`_PPP_LEGS`) ✓
- Backward-compat `NamedTuple` wrapper (phasemap/tests untouched) → Task 1 Steps 4–5 ✓
- Title as String or Observable → Task 1 ✓
- Row index + Z dict → Task 2 ✓
- Build-once network + slimmed callback → Task 3 ✓
- Tests: unchanged paths pass (Task 1 Step 5), plot-identity regression guard (Task 3), index unit tests (Task 2) → ✓

**Placeholder scan:** No TBD/TODO/"add error handling"/"similar to" — every code step shows full code. ✓

**Type consistency:** `_leg_drawspec` field names (`shaft`, `col`, `lw`, `tip`, `head`, `headrot`) are produced in Task 1 and consumed by the same task's `lift(d -> d[j].<field>, legdata)` calls — names match. `_PPP_LEGS` field names (`key`, `region`, `from`, `to`, `ctrl`) match their use in the `legdata` lift. `_build_row_index` key tuple `(atpase_frac, i_nadph, i_r5p)` matches the lookup in Task 3 (`row_index[(atp, in_, ir)]`, where `selected` is `(i_nadph, i_r5p)`). `_build_zdict` keyed by `atpase_frac` matches `zdict[atp]` in Task 3. ✓
