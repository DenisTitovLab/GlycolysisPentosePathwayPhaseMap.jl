# Resync GlycolysisPentosePathwayPhaseMap.jl to PentosePhosphatePathway.jl `main` — design

Date: 2026-10-13. Status: approved in brainstorming, pending written-spec review.

## Goal

The shipped grids were generated on 2026-06-18 (PhaseMap `0e92329`) against a PPP.jl that is now
47 commits behind `main` (`fcccae8`). Bring the package back in line with the current model and
make the metabolite panel compare like with like:

1. the recompute extension builds the model exactly as PPP.jl now defines it;
2. both shipped grids are regenerated from PPP.jl `main`;
3. the literature boxes in the metabolite panel sit on the model's pool basis
   (cytosolic water; free NADP(H)), with the conversion done in exactly one place;
4. every shipped data file records where it came from.

Success: `explore()` / `phasemap()` render the regenerated grids against the converted literature
table; the core test suite passes; `data/provenance.toml` names the PPP.jl and PPP_Experiments
commits behind the shipped data; docs and the local CLAUDE.md describe the new refresh workflow.

## Out of scope

- Marking the model's default operating point on the map; caching the R5Pase sink's ATP use.
- An old-vs-new grid comparison: the grids are regenerated and replaced.
- General registration (still gated on PPP.jl registration).

## §1 Environment and provenance

- `../PentosePhosphatePathway.jl` is on `main` at `fcccae8` and stays there (switched from
  `test/measured-keq`, which is pushed and untouched).
- Grids are regenerated in a scratch Julia environment that `Pkg.develop`s this package and
  `../PentosePhosphatePathway.jl`, so the solve uses exactly the checked-out commit.
- New tracked file `data/provenance.toml` with two tables, each owned by one writer that
  read-modify-writes only its own table (`TOML` stdlib; added to `[deps]`):
  - `[grids]` — written by `regenerate_grid`: `ppp_commit` (via `git -C pkgdir(PPP) rev-parse HEAD`
    when it is a git checkout, else `"v" * string(pkgversion(PPP))`), `ppp_dirty` (bool, git only),
    `cellmetabolism_version`, `generated` (date), the keyword arguments used, row counts per grid,
    and the count of cells whose `retcode != :Terminated` per grid.
  - `[literature]` — written by the PPP_Experiments export script (§4): `ppp_experiments_commit`,
    `ppp_experiments_dirty`, `source_files`, `nadp_seed`, `nadp_draws`, `generated`, row count.

## §2 Extension fix (`ext/PentosePhosphatePathwayExt.jl`)

- Initial state: `init0 = copy(PPP.glycolysis_ppp_init_conc)`. This replaces the hand-merge
  `merge_LArrays(PPP_init_conc, CellMetabolism.glycolysis_init_conc)`, which skips
  `PPP_glycolysis_init_overrides` (FIX-15/16/22: ADP, AMP, DHAP, Lactate, Pi, Pyruvate, NTP, NDP,
  F16BP, PEP) and mixes old and new adenylates.
- Parameters unchanged: `merge_LArrays(PPP.PPP_params, CellMetabolism.glycolysis_params)`
  (PPP.jl wins on shared keys).
- Supply ceilings, φ axes, and cached pool/flux columns unchanged (all names verified on `main`).
- `_regenerate_grid_impl` writes the `[grids]` table after both CSVs are written.

## §3 Grid regeneration

- Same axes as the shipped grids: 12 × 12 NADPH × R5P; 11 log-spaced ATP fractions 0.01–0.20;
  φ ∈ [1e-3, 1]; 2-D grid at `atpase_frac = 0.10` (= PPP.jl's default `0.2·HK1_cap`).
- Run in the background from the scratch environment; replaces `data/atp_demand_grid.csv` and
  `data/demand_grid.csv`.

## §4 Literature table

Conversion logic lives only in PPP_Experiments (`scripts/helper_files/metab_data.jl`, used by
`import_data.jl`); PhaseMap ships its output and does no conversion.

**PPP_Experiments** — new `scripts/phasemap_export/export_literature_pools.jl`:

- Includes `../helper_files/import_data.jl`, which already applies `apply_volume_basis!` (row
  `Basis` × species distribution class → cytosolic water) and the GAP/DHAP free-aldehyde
  corrections, and loads `free_pool_data`.
- **All rows are kept**, including drug-treated and glucose-pulse conditions
  (`<line>_<treatment>` Cell Types).
- Panel metabolites other than NADP/NADPH: the converted metabolomics column values
  (total-pool columns of `FREE_POOL_METABS` are skipped, as in `literature_pools`).
- NADPH: the cytosolic free-pool rows (`free_concentrations(free_pool_data, :NADPH)`).
- NADP (no direct free measurement): `nadp_draws` (default 1000) values of free NADPH / R, with
  NADPH resampled from the free NADPH rows and R from `draw_free_nadph_nadp_ratio` (the
  `docs/free_vs_total_measurements.md` §4 sampler), seeded (`nadp_seed`, default 20261013) so the
  file is reproducible. These rows carry `derived = true` and reference
  `"derived: free NADPH / R (free_vs_total_measurements.md §4)"`.
- Output `literature_pools.csv` (ASCII, long format), written to a path argument defaulting to
  `../GlycolysisPentosePathwayPhaseMap.jl/data/literature_pools.csv`; columns
  `metabolite, value_uM, reference, cell_type, derived`. Only the 16 panel metabolites are
  written. A metabolite with no measurements writes no rows and draws no box (PGLn today: the
  CSV has no PGLn values); the export errors only if a panel metabolite has no source at all
  (neither a metabolomics column nor a free-pool rule), which would signal a renamed column.
- Writes the `[literature]` table of the sibling `data/provenance.toml`.

**PhaseMap**:

- `metab_lit_path()` → `data/literature_pools.csv`; `data/Metabolite_Concentrations.csv` deleted.
- `load_metab_lit` reads the long format and returns the same `Dict{Symbol, Vector{Float64}}`
  (µM) keyed by `_PANEL_METAB_ORDER`. It is internal (not exported): no public API change.
- `_GAP_CORRECTION` / `_DHAP_CORRECTION` deleted.
- Panel y-label states the basis: `"[metabolite], µM (cytosolic water; NADP(H) free)"`.

## §5 Docs, version, CLAUDE.md

- `docs/src/recompute.md`: two-step refresh (grids via the extension; literature via the
  PPP_Experiments export) and `data/provenance.toml`.
- `docs/src/tutorial.md` / `science.md`: describe the literature boxes (model basis, free
  NADP(H), derived NADP, treated conditions included); recheck any statement about grid results
  against the regenerated grids.
- Re-render `docs/src/assets/explorer.png` with `scripts/render_readme_screenshot.jl`.
- `Project.toml` version 0.1.0 → 0.2.0 (shipped results change).
- Local `CLAUDE.md` (gitignored, not committed) rewritten to < 100 lines: general-mammalian
  framing, both refresh steps, provenance, hook-Ref extension pattern, backend import gotcha.

## §6 Testing

- Core suite (`Pkg.test()`):
  - `load_metab_lit` on the shipped table: returns every panel metabolite as a key; all values
    finite and > 0; every panel metabolite except PGLn non-empty; NADP and NADPH present.
  - The shipped CSV's NADP rows are all `derived = true` and no other metabolite's are.
  - `data/provenance.toml` parses; `[grids].ppp_commit` and `[literature].ppp_experiments_commit`
    are non-empty.
  - Both grids carry the expected pool/flux columns.
- Extension (manual; PPP.jl is unregistered so CI cannot load it): assert the extension's
  initial vector equals `PPP.glycolysis_ppp_init_conc`; run a 2 × 2 smoke `regenerate_grid` to
  scratch paths before the full run.

## §7 Git

- PhaseMap: branch `feat/ppp-main-resync` → PR to `main`.
- PPP_Experiments: branch `feat/phasemap-literature-export`; stage only this work's files.
