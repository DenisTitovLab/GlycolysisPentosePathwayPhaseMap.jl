# Split the explorer's metabolite panel into free ≈ total and free-pool axes — design

Date: 2026-10-14. Status: brainstormed with the user; all decisions below user-approved; one
review pass folded in. Pending written-spec review.

## Plain-English summary

The model tracks each metabolite as its **free** concentration: the molecules enzymes actually
see. Most literature measurements extract the cell and count the **total** pool, free plus
protein-bound. For most metabolites the two are about the same. For the cofactors ADP, AMP,
NAD⁺, NADH, NADP⁺ and NADPH most of the pool is bound, so only free-pool measurements are a fair
comparison.

PPP_Experiments' basic panel already splits its metabolite plot this way (spec
`PPP_Experiments/docs/superpowers/specs/2026-09-27-free-pool-metabolite-split-design.md`). This
change brings the same split to the `explore()` window: a left plot of free ≈ total metabolites
against extraction metabolomics, and a right plot of the six cofactors against free-pool data
only. Four cofactors (ADP, AMP, NAD⁺, NADH) are new to the explorer, so the shipped grids gain
four pool columns and are regenerated.

While designing this we found that the phase map's grids are solved with no F26BP, the
glycolysis activator, while the project's physiological operating point (and the basic panels)
use 5 µM. Without it, upper glycolysis runs near its limit and ADP builds up 20-fold. The
regenerated grids therefore also switch the phase map to 5 µM F26BP. PentosePhosphatePathway.jl's
own default stays at 0 for now.

## Goal

1. The explorer's metabolite row is two linked log-y axes: `_FREE_EQ_TOTAL_PANEL` (left) and
   `_FREE_POOL_PANEL` (right), with the same grouping as PPP_Experiments' `metab_data.jl`.
2. The right axis shows free-pool evidence as individual points in the basic panel's marker
   vocabulary; it never shows total-pool data.
3. Both shipped grids carry `ADP_uM, AMP_uM, NAD_uM, NADH_uM` and are solved at F26BP 5 µM.
4. The grouping cannot silently drift from PPP_Experiments: the export refuses to run on a
   mismatch.

Success: `explore()` opens on the physiological cell with both axes drawn; every right-axis
column has data points or reads "no data"; model values off the y-range show as labelled arrows;
the core test suite and the manual extension test pass; docs and the README screenshot show the
new panel; `[grids]` provenance records the F26BP value.

## Out of scope

- The basic panel's other free ≈ total metabolites (BPG, 3PG, 2PG, PEP, pyruvate, lactate, Pi):
  the grids do not store them. Adding them is the same kind of change, done separately.
- Free-ATP sensor rows: ATP stays on the left against metabolomics, as in the basic panel (its
  decision D1).
- The static `phasemap()` figure's layout (it has no metabolite panel). Its grids do change
  through D7.
- Changing PentosePhosphatePathway.jl's default F26BP (the user will do this separately).
- FIX-18 (adenylate starting values on the total basis). See "ADP and AMP" below: its effect is
  about 15% on every adenylate, not the 20-fold ADP gap.
- Any change to PPP_Experiments' basic panel or to `data/Free_Pool_Concentrations.csv`.

## Decisions

| # | Question | Decision |
|---|---|---|
| D1 | Right-axis scope | All six free-pool cofactors: ADP, AMP, NAD, NADH, NADP, NADPH. Grids are regenerated to carry the four new ones. |
| D2 | Right-axis data display | Individual points, basic-panel marker key: filled = measured; open = inferred, by us or by the original authors; circle = sensor-based; diamond = equilibrium-based. So ADP's direct equilibrium rows draw as open diamonds. |
| D3 | NADP literature values | Matched-pair derived points only (`derive_matched_free`), replacing the 1000-draw box. Today that is **one** point (see below). Confirmed by the user with the corrected count. |
| D4 | Where the grouping lives | A package constant, guarded by the PPP_Experiments export. Rejected: a `pool` column in the CSV (AMP has no rows to carry it); one shared plotting module for both repos (the two panels draw different things; only the grouping and marker key would be shared). |
| D5 | Left-axis display | Unchanged: static grey box-and-whisker per metabolite plus the purple model diamond. |
| D6 | Off-scale model values | Reactive labelled arrow at the floor or ceiling of the column, as in the basic panel, instead of silently dropping the diamond. Literature points get no off-scale handling; a test keeps them in range instead. |
| D7 | F26BP in the phase-map simulations | 5 µM (cytosolic water), the 2026-09-23 physiological operating point and the basic panels' value. Set in the extension, held in `_PHYSIOLOGICAL_DEMAND`, recorded in `[grids]`. PPP.jl's default (0) is untouched. |

**NADP (D3).** Free NADPH is measured only in HeLa (Tao 2017, Moon 2024), so only the HeLa
Sallin 2018 ratio (75.0) has a same-line partner. It gives 0.038 µM. The other three Sallin
lines and Veech 1969 are unmatched. This is what the basic panel's right axis shows. The old
1000-draw box spanned 1.2e-5 to 8.5 µM, nearly six decades, so it said little as a comparison
target. The tutorial states that the other ratio sources have no same-line free NADPH partner.

**ADP and AMP (D7).** The ADP/AMP split is set by F26BP, not by the adenylate starting values.
Physiological cell (NADPH φ 0.0432, R5P φ 0.1068, ATP 0.25 × 2·HK1_cap), same initial values,
solved 2026-10-14:

| F26BP | ATP | ADP | AMP | ATP + ADP + AMP |
|---|---|---|---|---|
| 0 (today's grids) | 6854 µM | 620 µM | 27.4 µM | 7501 µM |
| 5 µM (this change) | 7472 µM | 29.1 µM | 0.155 µM | 7501 µM |

The adenylate total is conserved, so the starting values only set the 7501 µM total. ADP's
total-basis starting value is 979 µM of it, so moving ADP and AMP to the free basis (FIX-18)
would scale every adenylate by roughly 15%. At 5 µM F26BP, model ADP sits inside the free-ADP
data (14–74 µM). No in-figure caveat is needed.

Free-pool data on the right axis today (`data/Free_Pool_Concentrations.csv`, checked
2026-10-14):

| Metabolite | Not derived (row's own method) | Derived (matched pair) |
|---|---|---|
| ADP | 6 equilibrium (open diamonds), 14–74 µM | — |
| AMP | — | — ("no data") |
| NAD | 6 sensor (filled circles), 44.6–106 µM | — |
| NADH | 1 sensor (filled circle), 0.13 µM | 1 equilibrium (open diamond), 0.34 µM |
| NADP | — | 1 sensor (open circle), 0.038 µM |
| NADPH | 2 sensor (filled circles), 2.6–3.1 µM | — |

## §1 Panel definition and operating point (`src/grid_io.jl`)

```julia
const _FREE_EQ_TOTAL_PANEL = (:Glucose, :G6P, :F6P, :F16BP, :DHAP, :GAP, :ATP,
                              :PGLn, :PGA, :Ru5P, :R5P, :X5P, :S7P, :E4P)
const _FREE_POOL_PANEL     = (:ADP, :AMP, :NAD, :NADH, :NADP, :NADPH)
const _PANEL_METAB_ORDER   = (_FREE_EQ_TOTAL_PANEL..., _FREE_POOL_PANEL...)
const _PHYSIOLOGICAL_DEMAND = (nadph_phi = 0.059, r5p_phi = 0.1068, atpase_frac = 0.25,
                               f26bp_uM = 5.0)
```

- The left order is today's order with NADP/NADPH removed. The right order is
  `FREE_POOL_METABS` in `metab_data.jl`.
- `_PANEL_METAB_ORDER` stays the single order used by `row_to_pools`, `load_metab_lit` and the
  export, so pools index as `[left; right]`. Its doc comment says it is the concatenation of the
  two panels.
- The comment on the two panel tuples names `PPP_Experiments/scripts/helper_files/metab_data.jl`
  as the grouping's source and the export guard (§3) as its enforcement.
- The `_PHYSIOLOGICAL_DEMAND` comment gains the F26BP source (PPP_Experiments
  `notes/2026-09-23_metabolite-dynamics-physiological-operating-point.md`) and states that it
  differs from PPP.jl's default of 0.

## §2 Grids (`ext/PentosePhosphatePathwayExt.jl`, `data/*.csv`, `src/grid_io.jl`)

Extension:
- `_DEMAND_POOLS` appends `:ADP, :AMP, :NAD, :NADH`. These are model species
  (`glycolysis_ppp_init_conc` carries them) and are read from `sol.u[end]` exactly as the
  existing pools are (× 1e6 → µM; NaN on solver-error rows).
- New keyword `f26bp_uM = GPPM._PHYSIOLOGICAL_DEMAND.f26bp_uM`, threaded through
  `_regenerate_grid_impl` → `run_atp_grid` / `run_demand_grid` → `_demand_setup`, which sets
  `init0.F26BP = f26bp_uM * 1e-6` after the copy. The value is cytosolic water, so there is no
  `whole_cell_to_cytosolic` conversion, matching `Exp_basic_panel_v2.jl:66`. F26BP has no
  reactions, so it stays at this value.
- `_grid_provenance` records `f26bp_uM`. A `[grids]` table without the key means 0 (today's
  grids).

Regeneration runs in the background from the PPP_Experiments environment, against PPP.jl `main`
(`fcccae8`, the commit behind today's grids), with the other keyword arguments as recorded in
`[grids]`. There are two runs, both written to scratch:

1. **Drift check, `f26bp_uM = 0`.** Compare with today's grids, column by column by name (the
   column order changes because the new pools are inserted). Use `isequal` for numbers (exact,
   NaN-aware) and `==` for the `retcode` and `mode` strings, and check that the row order
   matches. The expected result is bit-identical, since the solve is deterministic and
   unchanged. Any difference means PPP.jl or the environment moved: stop and report. This run
   isolates the new columns from the F26BP change and is not shipped.
2. **Ship run, `f26bp_uM = 5.0`.** Summarise the change from today's grids: collapsed and
   non-converged counts per ATP level, how many cells change `mode`, and the cycle-index shift
   at the physiological cell. Also re-check the reviewer's note that the lowest-NADPH-demand
   cells are ATP-collapsed at ATP 0.25 under F26BP 0. The summary goes in the package commit
   body and the PR description.

Then copy the ship-run CSVs into `data/` and copy its `[grids]` table into
`data/provenance.toml` with `GPPM.update_provenance!`, which leaves `[literature]` untouched.

Viewer-side column checks:
- `load_grid` checks only the columns every renderer needs: `retcode`, `cycle_index`, `mode`,
  `i_nadph`, `i_r5p`, `nadph_phi`, `r5p_phi`, and `"$(k)_uM_min"` for every flux key the
  renderers read. A missing column throws an `ArgumentError` that names it.
- `build_explorer` also checks `atpase_frac` and `"$(s)_uM"` for `_PANEL_METAB_ORDER`. A missing
  pool column throws "grid predates the cofactor pool columns (…); regenerate it". `phasemap()`
  never reads pools, so it keeps working on an older grid.

## §3 Literature table (PPP_Experiments `scripts/phasemap_export/`)

The schema of `data/literature_pools.csv` becomes
`metabolite, value_uM, reference, cell_type, derived, method`, where `method` ∈
`metabolomics`, `sensor`, `equilibrium`.

`literature_pools_table(metab_df, free_df, free_derived, panel)`:
- `FREE_EQ_TOTAL_METABS` entries: every metabolomics row, as today, with `derived = false` and
  `method = "metabolomics"`. A panel entry with no metabolomics column still throws "no
  literature source", as today.
- `FREE_POOL_METABS` entries: first the direct cytosolic concentration rows (`derived = false`,
  the row's `Method`, `Value × 1e6`). Then the `free_derived` rows for that metabolite
  (`derived = true`, the ratio's `Method`, `Value × 1e6`, the matched `Cell Type`). Their
  `reference` names both inputs:
  `"derived: <ratio Reference> x median free <partner> (n=<n_partner_rows>)"`. `free_derived` is
  the `free_pool_derived` global that `import_data.jl` already builds with
  `derive_matched_free`; the export does not re-derive. A free-pool metabolite with no rows
  (AMP) is not an error.
- Removed: the NADP 1000-draw branch; `NADP_DRAWS`, `NADP_SEED`, `NADP_DERIVED_REFERENCE`,
  `NADP_DERIVED_CELL_TYPE`; the `nadp_draws`/`nadp_seed` keywords; `using Random` in the export
  and its test; the "direct free NADP rows exist" guard (direct rows are now simply included);
  and the `nadp_seed`/`nadp_draws` keys of `[literature]`. The §4 ratio sampler stays: it is
  still used by `bootstrap_metab_data`.
- Kept: the finite-and-positive check and the ASCII check.

**Grouping guard.** `check_panel_grouping(free_eq_total_panel, free_pool_panel)` lives in
`literature_export.jl`. It throws unless `free_pool_panel == FREE_POOL_METABS` (same order) and
every `free_eq_total_panel` entry is in `FREE_EQ_TOTAL_METABS`.
`export_literature_pools.jl` calls it with `GPPM._FREE_EQ_TOTAL_PANEL` and
`GPPM._FREE_POOL_PANEL` before anything is written.

## §4 Viewer (`src/metab_panel.jl`, `src/explorer.jl`)

**Loader.** `load_metab_lit(path; order)` returns
`Dict{Symbol, Vector{@NamedTuple{value::Float64, derived::Bool, method::Symbol}}}`, with a
vector (possibly empty) for every `order` entry. It throws on a `method` outside the three
allowed values, and on a left-panel metabolite with a non-`metabolomics` row or a right-panel
metabolite with a `metabolomics` row.

**Marker style.** `_free_pool_point_style(derived, method)` is copied from
`PPP_Experiments/scripts/helper_files/plot_metab_split.jl` (`free_pool_point_style`), with a
comment naming the source. Both panels then use one marker key.

**`draw_metab_split!(gl, lit, metab_obs) -> (; ax_left, ax_right)`** replaces
`draw_metab_panel!`. `gl` is a `GridLayout`; `metab_obs` is the length-20
`_PANEL_METAB_ORDER` Observable, split as `[1:14]` left and `[15:20]` right.

Layout:
- Two `Axis` in `gl[1, 1]` and `gl[1, 2]`, plus a `Legend` in `gl[1, 3]`, outside both data
  areas. The basic-panel mockup showed an in-axis legend hiding the NADP data.
- Column widths `Auto(14)` and `Auto(6)`, `colgap!` 12. `linkyaxes!`. Log10 y, with
  `ylims!(1e-4, 1e5)` set before switching to log10, as today. The right axis hides its y tick
  labels but keeps ticks and grid.
- **x-limits pinned** on both axes: `xlims!(ax, 0.5, n + 0.8)`. The extra 0.3 on the right
  leaves room for an off-scale label in the last column. The axes never auto-limit x, because
  the right axis's marks are sparse and reactive; unpinned, the review render clipped the NADPH
  column.
- Titles: left "free ≈ total (extraction metabolomics)", right "free pool (free-pool data
  only)". Check the screenshot for overflow. One y-label: "[metabolite], µM (cytosolic water)".
- A right-axis column with no data gets a second tick-label line, "AMP\nno data".
- Why the floor drops from 1e-3 to 1e-4: in today's grids, NADP falls to 3.6e-4 µM and R5P to
  3.1e-4 µM in cells that have not collapsed. With 1e-4, no live cell has any pool off-scale;
  only ATP-collapsed cells draw arrows, up to about 13 per cell. Re-check these numbers on the
  F26BP 5 µM grids and update the sentence in the docstring.

Marks:
- Left: today's static grey boxes from the `value`s, plus purple diamonds.
- Right: static points, one per row, spread over ±0.2 of the column when there are several,
  styled by `_free_pool_point_style`.
- Both axes: purple model diamonds lifted from `metab_obs` for in-range values. Each axis also
  has one lifted off-scale `scatter!` (`:dtriangle` at the floor, `:utriangle` at the ceiling,
  purple) and one lifted `text!` giving the true value (`%.1e`). Zero and tiny negative solver
  values count as below the floor, so a vanished AMP pool stays visible. Only non-finite values
  (solver-error rows) draw nothing. The plot count per axis stays fixed when the selection
  changes; the review confirmed this works for a varying off-scale count.
- Legend: "literature (metabolomics)" box, "selected cell" diamond, "selected cell, off scale"
  `:dtriangle`, and the three free-pool marker kinds.

`explorer.jl`: `fig[4, 1:3]` holds a `GridLayout` passed to `draw_metab_split!`, and
`rowsize!(fig.layout, 4, Relative(0.30))` stays. Zoom and pan interactions are deregistered on
both axes. Changing the selection only updates Observables. `_METAB_YLABEL` is kept (tests find
the left axis by it), with the new text.

## §5 Tests

**Package, removed or replaced** (fail after this change):
- `test_metab_panel.jl:12`: every non-PGLn panel metabolite has data. Replace with: every panel
  metabolite except PGLn and AMP has data.
- `test_metab_panel.jl:14`: the column names. Add `method`.
- `test_metab_panel.jl:16`: no `derived` rows outside NADP. Replace with: `derived` only on
  right-panel metabolites, and all NADP rows are `derived`.
- `test_metab_panel.jl:33-44`: the `draw_metab_panel!` testset, including the "free" check on
  `_METAB_YLABEL`. Replace with the split testset below.
- `test_metab_panel.jl:5-31`: both `load_metab_lit` testsets. Update for the NamedTuple return
  and the `method` validation.

**Package, new:**
- `test_grid_io.jl`:
  - the shipped grids carry all 20 `<sym>_uM` columns;
  - the panel tuples partition `_PANEL_METAB_ORDER`, with no overlap;
  - `load_grid` throws on a missing flux column;
  - `build_explorer` throws a clear error on a grid missing `ADP_uM`;
  - `phasemap` still renders a grid that lacks the pool columns;
  - `[grids]` provenance has `f26bp_uM == _PHYSIOLOGICAL_DEMAND.f26bp_uM`.
- `test_metab_panel.jl`:
  - `method` takes only the three values;
  - left-panel rows are all `metabolomics` and right-panel rows never are;
  - AMP has no rows;
  - every shipped value lies inside the y-limits (this guards D6's "no off-scale data");
  - `draw_metab_split!` renders offscreen with a synthetic `lit` that includes an off-scale
    value, a zero, a NaN and an empty column;
  - it re-renders after Observable updates with the same plot count on each axis;
  - `ax.finallimits[]` x-range equals the pinned range after an update.
- `test_explorer.jl`: finds both metabolite axes, with interactions deregistered on both.
- `test/manual/test_extension.jl`: re-pin the existing numeric expectations
  (`30 < V_R5Pase < 40`, ATP > 1000, etc.) to F26BP 5 µM values. Require finite, positive
  ADP/AMP/NAD/NADH. Check that `f26bp_uM = 0` reproduces a cell of today's grids.

**PPP_Experiments** (`scripts/phasemap_export/test_literature_export.jl`): replace the
200-draw, seed and direct-NADP-guard tests (lines 34-56, 70-72). New tests cover:
- the new schema and `method` per group;
- direct plus derived rows for a free-pool metabolite, with the two-input derived `reference`;
- no error for a free-pool metabolite with no rows;
- `check_panel_grouping` rejecting a reordered or mismatched tuple.

## §6 Docs and versioning

- `docs/src/tutorial.md`, "The metabolite panel":
  - rewrite for two axes, 20 metabolites, the marker key and the "no data" column;
  - say NADP is one matched point and why;
  - drop the derived-NADP-box text.
- `docs/src/tutorial.md` and `science.md`: state the operating point's F26BP (5 µM) where the
  physiological demand is described.
- `docs/src/recompute.md`:
  - "free-pool NADP(H)" becomes "free-pool cofactors";
  - document the `method` column;
  - drop the seed and draw count;
  - document `f26bp_uM`.
- README and `docs/src/assets/explorer.png`: regenerate with `scripts/render_readme_screenshot.jl`.
- Package `CLAUDE.md` is gitignored, so this is a local edit only: note the two panel tuples,
  the export guard and `f26bp_uM`, and remove the stale `_GAP_CORRECTION`/`_DHAP_CORRECTION`
  text.
- Version 0.2.0 → 0.3.0: the grid and literature CSV schemas change, the shipped grids change
  (F26BP), and `explore(grid)` rejects grids without the cofactor pools.

## §7 Order of work and repos

Each repo gets one feature branch, squash-merged. The plan asks the user before creating either
branch. Stage only this work's files: PPP_Experiments has unrelated uncommitted work.

1. **Package branch:**
   - add the panel tuples, `f26bp_uM`, `_DEMAND_POOLS` and the column checks;
   - run both regenerations (§2);
   - install the ship-run grids.
2. **PPP_Experiments branch:**
   - make the export changes and the guard (§3) against the package branch's tuples (the
     environment `develop`s the package by path);
   - merge that PR.
3. **From PPP_Experiments `main` at the squash commit:** run the export into the package
   branch's `data/`. `[literature]` then records a clean commit that will not be orphaned, which
   avoids the follow-up fix the resync needed (`458e79f`).
4. **Package branch:**
   - the viewer split (§4), tests (§5), docs and screenshot (§6);
   - merge.

The package suite can be red between steps 1 and 3: the tuples list metabolites that the
literature table does not have yet. Only the branch tip has to be green, since it is
squash-merged.
