# Free-Pool Split of the Explorer's Metabolite Panel — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Split the `explore()` metabolite panel into a free ≈ total axis (metabolomics boxes) and a free-pool axis (ADP, AMP, NAD, NADH, NADP, NADPH against free-pool points). Regenerate the grids with the four new cofactor pools at F26BP 5 µM.

**Architecture:** Two repos.
- **PhaseMap** (`/home/james/projects/GlycolysisPentosePathwayPhaseMap.jl`, branch `feat/free-pool-panel-split`) gains:
  - an F26BP keyword in the recompute extension;
  - four new pool columns;
  - two panel tuples;
  - grid column checks;
  - a two-axis Makie panel.
- **PPP_Experiments** (`/home/james/projects/PPP_Experiments`, branch `feat/phasemap-free-pool-export`) changes the literature export: a `method` column, matched-pair rows instead of the 1000 NADP draws, and a guard that the two repos' grouping agrees.

Every model run and every export uses the PPP_Experiments environment, which `develop`s PhaseMap and `../PentosePhosphatePathway.jl` by path.

**Tech Stack:** Julia 1.12 (package compat 1.10), Makie 0.24 / CairoMakie / GLMakie, CSV.jl, DataFrames.jl, TOML stdlib, OrdinaryDiffEq `Rodas5P`.

**Spec:** `docs/superpowers/specs/2026-10-14-free-pool-panel-split-design.md` (PhaseMap repo, committed `a3a70f6`). Read it before starting any task.

## Global Constraints

- **Repos and state:**
  - PPP.jl checkout `/home/james/projects/PentosePhosphatePathway.jl` must be on `main` at `fcccae8a196ce9cfba23fae91ca77c4be4540807`, clean, before any model run. Never switch it.
  - Solver settings are unchanged: `Rodas5P()`, `TerminateSteadyState(1e-13, 1e-8)`, `abstol = 1e-14`, `reltol = 1e-9`.
  - Grid axes are unchanged (they are the extension's defaults): `n_nadph = 12`, `n_r5p = 12`, `atp_fracs = 10 .^ range(log10(0.01), log10(0.25), length = 11)`, `atpase_frac = 0.25`, `phi_lo = 1e-3`, `phi_hi = 1.0`.
- **Panel and operating point:**
  - `_FREE_EQ_TOTAL_PANEL = (:Glucose, :G6P, :F6P, :F16BP, :DHAP, :GAP, :ATP, :PGLn, :PGA, :Ru5P, :R5P, :X5P, :S7P, :E4P)`.
  - `_FREE_POOL_PANEL = (:ADP, :AMP, :NAD, :NADH, :NADP, :NADPH)`.
  - `_PANEL_METAB_ORDER = (_FREE_EQ_TOTAL_PANEL..., _FREE_POOL_PANEL...)`.
  - `_PHYSIOLOGICAL_DEMAND = (nadph_phi = 0.059, r5p_phi = 0.1068, atpase_frac = 0.25, f26bp_uM = 5.0)`.
  - F26BP is set as `init0.F26BP = f26bp_uM * 1e-6`, with no volume conversion (cytosolic water).
- **Data files:**
  - `literature_pools.csv` columns are exactly `metabolite, value_uM, reference, cell_type, derived, method`, with `method ∈ {metabolomics, sensor, equilibrium}`. ASCII only.
  - The metabolite y-limits are `(1e-4, 1e5)` µM. The y-label is exactly `"[metabolite], µM (cytosolic water)"`.
  - `data/provenance.toml`: `[grids]` is written only from the extension's provenance, and `[literature]` only by the PPP_Experiments export.
- **Versioning:** PhaseMap version 0.2.0 → 0.3.0.
- **Naming and style:**
  - 6-phosphogluconate is `PGA` everywhere; never `6PG`.
  - Section banners are plain ASCII `#`, 90 characters wide, with the title centred.
- **Git:**
  - Stage named files only; never `git add -A`. Run `git branch --show-current` before every commit.
  - Commit subjects use `type(scope): summary`. A body that a non-coder couldn't follow opens with a plain-English summary.
  - Every commit message ends with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Environment:** GLMakie needs an X display even to load. If a command that loads PhaseMap fails with a GLFW or display error, prefix it with `xvfb-run -a`.
- **Local files:** PhaseMap's `CLAUDE.md` is gitignored. Edit it, never commit it.
- **Scratch:** use `SCRATCH=/tmp/claude-1000/-home-james-projects-PPP-Experiments/a3d4eca5-2f02-4011-8c30-18a62260a8ba/scratchpad/regen`.
- **Stop rule:** never loosen a numeric test bound to make it pass. If a pinned model number moves, stop and report the old and new values.

## Review Focus

1. **An older grid opened with `explore(grid = old.csv)`**, i.e. one without the cofactor pool columns. Expected: an `ArgumentError` that names the missing columns and says to regenerate, never a crash mid-draw. Meanwhile `phasemap(grid = old.csv)` still renders. Pinned in Task 3 (both tests).
2. **A selected cell with collapsed ATP**, where many pools fall below 1e-4 µM or reach exactly 0 (AMP). Expected: a labelled purple triangle at the floor for each such pool. Nothing vanishes and nothing throws on log-y. Pinned in Task 6 by `_offscale_marks` unit tests (zero, negative, NaN, above the ceiling) and the render test.
3. **Reselecting cells changes how many marks are off-scale** (3 → 0 → 5). Expected: no error, and a fixed plot count per axis. Pinned in Task 6's Observable-update test.
4. **The two repos' grouping drifts apart**, for example someone adds `:Pi` to `FREE_POOL_METABS`. Expected: the export refuses to run and names the mismatch. Pinned in Task 4's `check_panel_grouping` tests.
5. **A literature export that puts a total-pool row on the free-pool axis, or a value outside the y-limits.** Expected: `load_metab_lit` throws (wrong side), and the shipped-table test fails (out of range). Pinned in Task 5 (shipped table) and Task 6 (loader validation).

---

## File Structure

PhaseMap:
- `src/grid_io.jl`: the panel tuples, `_PHYSIOLOGICAL_DEMAND.f26bp_uM`, the column-check helper and the `load_grid` check (Tasks 1, 3).
- `ext/PentosePhosphatePathwayExt.jl`: `_DEMAND_POOLS`, the `f26bp_uM` keyword and its provenance key (Task 1).
- `data/atp_demand_grid.csv`, `data/demand_grid.csv`, `data/provenance.toml`: regenerated (Task 2); `[literature]` re-exported (Task 5).
- `data/literature_pools.csv`: re-exported (Task 5).
- `src/metab_panel.jl`: the loader, marker style, off-scale marks and `draw_metab_split!`, replacing `draw_metab_panel!` (Tasks 5, 6).
- `src/explorer.jl`: the pool-column check and the nested `GridLayout` wiring (Tasks 3, 6).
- `test/test_grid_io.jl`, `test/test_metab_panel.jl`, `test/test_explorer.jl`, `test/test_phasemap.jl`, `test/manual/test_extension.jl`.
- `docs/src/tutorial.md`, `docs/src/science.md`, `docs/src/recompute.md`, `docs/src/assets/explorer.png`, `Project.toml`, and the local `CLAUDE.md` (Task 7).

PPP_Experiments:
- `scripts/helper_files/metab_data.jl`: `derive_matched_free` gains `Ratio` and `Partner` columns (Task 4).
- `scripts/helper_files/test/test_free_pool.jl`: tests for the new columns (Task 4).
- `scripts/phasemap_export/literature_export.jl`, `export_literature_pools.jl`, `test_literature_export.jl` (Task 4).

---

### Task 1: Extension — F26BP keyword and four new pool columns

**Repo:** PhaseMap, branch `feat/free-pool-panel-split`.

**Files:**
- Modify: `src/grid_io.jl:9-15` (`_PHYSIOLOGICAL_DEMAND` and its comment)
- Modify: `ext/PentosePhosphatePathwayExt.jl:28-30, 44-57, 106-158, 163-182`
- Test: `test/manual/test_extension.jl`

**Interfaces:**
- Produces: `GPPM._PHYSIOLOGICAL_DEMAND.f26bp_uM == 5.0`.
- Produces: `ext._demand_setup(atpase_frac, n_nadph, n_r5p, phi_lo, phi_hi; f26bp_uM = GPPM._PHYSIOLOGICAL_DEMAND.f26bp_uM)`.
- Produces: the keyword `f26bp_uM` on `run_demand_grid`, `run_atp_grid` and `regenerate_grid` / `_regenerate_grid_impl`.
- Produces: the `[grids]` key `"f26bp_uM"::Float64`.
- Produces: grid CSV columns `ADP_uM, AMP_uM, NAD_uM, NADH_uM`, appended after `ATP_uM`.

- [ ] **Step 1: Write the failing manual tests**

In `test/manual/test_extension.jl`, replace the three `init0` assertions in the `"recompute extension"` testset:

```julia
    # initial state = PPP.jl's canonical combined init, except F26BP at the phase map's
    # physiological 5 µM (cytosolic water; PPP.jl's own default is 0)
    s = ext._demand_setup(0.10, 2, 2, 1e-3, 1.0)
    @test keys(s.init0) == keys(PPP.glycolysis_ppp_init_conc)
    @test s.init0.F26BP == G._PHYSIOLOGICAL_DEMAND.f26bp_uM * 1e-6 == 5e-6
    same = [k for k in keys(s.init0) if k != :F26BP]
    @test all(getproperty(s.init0, k) == getproperty(PPP.glycolysis_ppp_init_conc, k) for k in same)
    @test s.init0 !== PPP.glycolysis_ppp_init_conc          # a copy, never the package global
    @test PPP.glycolysis_ppp_init_conc.F26BP == 0.0         # the package global is untouched
    s0 = ext._demand_setup(0.10, 2, 2, 1e-3, 1.0; f26bp_uM = 0.0)
    @test collect(s0.init0) == collect(PPP.glycolysis_ppp_init_conc)   # 0 reproduces PPP.jl
```

In the same testset, after `@test prov["atpase_frac"] == ...`, add:

```julia
    @test prov["f26bp_uM"] == G._PHYSIOLOGICAL_DEMAND.f26bp_uM
    for s in (:ADP, :AMP, :NAD, :NADH)
        @test "$(s)_uM" in names(atp)
    end
```

At the end of the `"cold start, flux units, physiological defaults"` testset, add:

```julia
    # the four cofactor pools are cached, finite and positive at the physiological cell
    # (AMP included: a vanished AMP is what the panel's below-floor arrow is for)
    for k in (:ADP_uM, :AMP_uM, :NAD_uM, :NADH_uM)
        v = getproperty(r, k)
        @test isfinite(v) && v > 0
    end
    # F26BP sets the adenylate split (spec D7): 5 µM keeps ADP in the free-ADP data range
    @test 10 < r.ADP_uM < 100
```

- [ ] **Step 2: Run the manual test to verify it fails**

Run from the PhaseMap root:

```bash
git -C /home/james/projects/PentosePhosphatePathway.jl rev-parse HEAD   # must print fcccae8a196ce9cfba23fae91ca77c4be4540807
git -C /home/james/projects/PentosePhosphatePathway.jl status --porcelain  # must print nothing
julia --project=/home/james/projects/PPP_Experiments test/manual/test_extension.jl
```

Expected: FAIL. `_PHYSIOLOGICAL_DEMAND` has no field `f26bp_uM`, and `_demand_setup` has no keyword `f26bp_uM`.

- [ ] **Step 3: Add `f26bp_uM` to the operating point**

In `src/grid_io.jl`, replace the `_PHYSIOLOGICAL_DEMAND` comment and constant (lines 9-15):

```julia
# Default operating point: PPP.jl's default sink settings (src/model_parameters.jl:
# NADPHox_Vmax = 0.059·(6·RPI_cap), R5Pase_Vmax = 0.1068·(1.5·RPI_cap)) as fractions φ of the
# same supply ceilings the grids use, at an ATPase demand of 0.5·HK1_cap (= 0.25 of the
# 2·HK1_cap ATP ceiling; PPP.jl's own default is 0.2·HK1_cap). F26BP is 5 µM (cytosolic water),
# the physiological operating point of PPP_Experiments
# notes/2026-09-23_metabolite-dynamics-physiological-operating-point.md and the value its basic
# panels use; PPP.jl's own default is 0, at which upper glycolysis runs near its limit and ADP
# builds up ~20-fold. The 2-D grid is solved at this ATP level, and phasemap()/explore() open
# on the grid cell nearest it. The manual extension test pins the two φ against PPP.jl.
const _PHYSIOLOGICAL_DEMAND = (nadph_phi = 0.059, r5p_phi = 0.1068, atpase_frac = 0.25,
                               f26bp_uM = 5.0)
```

- [ ] **Step 4: Thread `f26bp_uM` through the extension and add the pools**

In `ext/PentosePhosphatePathwayExt.jl`:

Replace `_DEMAND_POOLS` (lines 29-30):

```julia
const _DEMAND_POOLS = (:R5P, :S7P, :NADPH, :NADP, :Ru5P, :X5P, :E4P, :F6P, :GAP, :G6P,
                       :Glucose, :F16BP, :DHAP, :PGLn, :PGA, :ATP, :ADP, :AMP, :NAD, :NADH)
```

Replace `_demand_setup` (lines 43-57):

```julia
# Shared per-sweep state: base params/init, the NADPH/R5P supply ceilings, the demand-φ axes.
function _demand_setup(atpase_frac, n_nadph, n_r5p, phi_lo, phi_hi;
                       f26bp_uM = GPPM._PHYSIOLOGICAL_DEMAND.f26bp_uM)
    params0 = PPP.merge_LArrays(PPP.PPP_params, CellMetabolism.glycolysis_params)
    # PPP.jl's canonical combined init: PPP_init_conc, then PPP_glycolysis_init_overrides
    # (FIX-15/16/22), then CellMetabolism's. A hand-merge skips the overrides. Copied: the
    # export is a package-level array that the sweep must never write into.
    init0   = copy(PPP.glycolysis_ppp_init_conc)
    # F26BP is cytosolic water (no whole_cell_to_cytosolic), as in PPP_Experiments'
    # Exp_basic_panel_v2.jl. It has no reactions, so it holds this value at steady state.
    init0.F26BP = f26bp_uM * 1e-6
    RPI_cap = params0.RPI_Vmax * params0.RPI_Conc
    S_NADPH = 6 * RPI_cap                  # 2·V_ox,max, V_ox,max = 3·RPI_cap (Ru5P-disposal limit)
    S_R5P   = 1.5 * RPI_cap                # RPI fwd + TKT_Rxn1-reverse bypass; 2/3 of R5P via RPI
    params0.ATPase_Vmax = atpase_frac * (params0.HK1_Vmax * params0.HK1_Conc * 2)  # φ·S_ATP
    nadph_phis = 10 .^ range(log10(phi_lo), log10(phi_hi), length = n_nadph)
    r5p_phis   = 10 .^ range(log10(phi_lo), log10(phi_hi), length = n_r5p)
    return (; params0, init0, S_NADPH, S_R5P, nadph_phis, r5p_phis)
end
```

Replace `run_demand_grid` and `run_atp_grid` (lines 106-137):

```julia
# 2-D NADPH × R5P demand grid at a fixed ATP background.
function run_demand_grid(; n_nadph = 12, n_r5p = 12,
                           atpase_frac = GPPM._PHYSIOLOGICAL_DEMAND.atpase_frac,
                           phi_lo = 1e-3, phi_hi = 1.0,
                           f26bp_uM = GPPM._PHYSIOLOGICAL_DEMAND.f26bp_uM)
    setup = _demand_setup(atpase_frac, n_nadph, n_r5p, phi_lo, phi_hi; f26bp_uM)
    rows = NamedTuple[]
    for (ir, rphi) in enumerate(setup.r5p_phis)
        append!(rows, _demand_column(setup, ir, rphi))
    end
    return DataFrame(rows)
end

# ATP-demand levels of the 3-D grid: 11 log-spaced fractions of the 2·HK1_cap ATP ceiling from
# 0.01 up to the physiological operating point's 0.25 (= 0.5·HK1_cap).
_default_atp_fracs() =
    10 .^ range(log10(0.01), log10(GPPM._PHYSIOLOGICAL_DEMAND.atpase_frac), length = 11)

# 3-D ATP × NADPH × R5P grid: the demand grid repeated across log-spaced ATP-demand levels,
# each block tagged with its `atpase_frac`. Serial (run-once, offline) — minutes for a full grid.
function run_atp_grid(; n_nadph = 12, n_r5p = 12, atp_fracs = _default_atp_fracs(),
                        phi_lo = 1e-3, phi_hi = 1.0,
                        f26bp_uM = GPPM._PHYSIOLOGICAL_DEMAND.f26bp_uM)
    frames = DataFrame[]
    for f in atp_fracs
        s = _demand_setup(f, n_nadph, n_r5p, phi_lo, phi_hi; f26bp_uM)
        for ir in 1:n_r5p
            sub = DataFrame(_demand_column(s, ir, s.r5p_phis[ir]))
            sub.atpase_frac .= f
            push!(frames, sub)
        end
    end
    return vcat(frames...)
end
```

In `_regenerate_grid_impl` (lines 142-158), add the keyword `f26bp_uM = GPPM._PHYSIOLOGICAL_DEMAND.f26bp_uM` after `phi_hi = 1.0`. Pass `f26bp_uM` to `run_atp_grid`, `run_demand_grid` and `_grid_provenance`:

```julia
function _regenerate_grid_impl(; atp_out = GPPM.default_atp_grid(),
                                 demand_out = GPPM.default_demand_grid(),
                                 provenance_out = joinpath(dirname(atp_out), "provenance.toml"),
                                 n_nadph = 12, n_r5p = 12,
                                 atp_fracs = _default_atp_fracs(),
                                 atpase_frac = GPPM._PHYSIOLOGICAL_DEMAND.atpase_frac,
                                 phi_lo = 1e-3, phi_hi = 1.0,
                                 f26bp_uM = GPPM._PHYSIOLOGICAL_DEMAND.f26bp_uM)
    mkpath(dirname(atp_out)); mkpath(dirname(demand_out))
    atp = run_atp_grid(; n_nadph, n_r5p, atp_fracs, phi_lo, phi_hi, f26bp_uM)
    CSV.write(atp_out, atp)
    demand = run_demand_grid(; n_nadph, n_r5p, atpase_frac, phi_lo, phi_hi, f26bp_uM)
    CSV.write(demand_out, demand)
    GPPM.update_provenance!(provenance_out, "grids",
        _grid_provenance(atp, demand; n_nadph, n_r5p, atp_fracs, atpase_frac, phi_lo, phi_hi,
                         f26bp_uM))
    @info "Regenerated grids" atp_out demand_out provenance_out atp_rows=nrow(atp) demand_rows=nrow(demand)
    return (; atp_out, demand_out, provenance_out)
end
```

In `_grid_provenance` (line 163), add `f26bp_uM` to the keyword list. In the `Dict`, after `"phi_lo" => ..., "phi_hi" => ...`, add `"f26bp_uM" => Float64(f26bp_uM),`.

In `src/GlycolysisPentosePathwayPhaseMap.jl`, extend the `regenerate_grid` docstring with one sentence: "`f26bp_uM` (default `_PHYSIOLOGICAL_DEMAND.f26bp_uM`, 5 µM) sets the F26BP initial and steady-state value, in µM of cytosolic water; `[grids]` records it."

- [ ] **Step 5: Run the manual test to verify it passes**

Run: `julia --project=/home/james/projects/PPP_Experiments test/manual/test_extension.jl`

Expected: PASS for every testset. If `30 < r.V_R5Pase_uM_min < 40`, `col[7].ATP_uM > 1000` or `10 < r.ADP_uM < 100` fails, **stop**. Report the value; do not change the bound.

- [ ] **Step 6: Run the core suite (unchanged grids, still green)**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`

Expected: PASS. Grids and panel are unchanged so far.

- [ ] **Step 7: Commit**

```bash
git branch --show-current   # feat/free-pool-panel-split
git add src/grid_io.jl src/GlycolysisPentosePathwayPhaseMap.jl ext/PentosePhosphatePathwayExt.jl test/manual/test_extension.jl
git commit -F - <<'EOF'
feat(ext): solve the grids at 5 uM F26BP and cache ADP, AMP, NAD and NADH

Plain English: the phase-map simulations now run with the glycolysis activator F26BP at its
physiological 5 uM, as the PPP_Experiments basic panels do, instead of PentosePhosphatePathway.jl's
default of zero. The grids will also record four more cofactor concentrations. The shipped grids
are regenerated in the next commit.

f26bp_uM keyword on _demand_setup / run_*_grid / regenerate_grid (default
_PHYSIOLOGICAL_DEMAND.f26bp_uM = 5.0), recorded in [grids]; f26bp_uM = 0 reproduces PPP.jl's
init exactly. _DEMAND_POOLS gains ADP, AMP, NAD, NADH.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 2: Regenerate the grids (drift check at F26BP 0, ship run at 5 µM)

**Repo:** PhaseMap, branch `feat/free-pool-panel-split`.

**Files:**
- Replace: `data/atp_demand_grid.csv`, `data/demand_grid.csv`
- Modify: `data/provenance.toml` (`[grids]` table only)
- Test: `test/test_grid_io.jl`

**Interfaces:**
- Consumes: `regenerate_grid(; atp_out, demand_out, f26bp_uM)` from Task 1.
- Produces: shipped grids with `ADP_uM, AMP_uM, NAD_uM, NADH_uM`, solved at F26BP 5 µM; `[grids]["f26bp_uM"] == 5.0`.

- [ ] **Step 1: Write the failing tests**

In `test/test_grid_io.jl`, add this testset after `"shipped grids carry the viewer's columns"`:

```julia
@testset "shipped grids carry the cofactor pools and the F26BP setting" begin
    for path in (G.default_atp_grid(), G.default_demand_grid())
        cols = names(CSV.read(path, DataFrame; limit = 1))
        for s in (:ADP, :AMP, :NAD, :NADH)
            @test "$(s)_uM" in cols
        end
    end
    @test G.read_provenance()["grids"]["f26bp_uM"] == G._PHYSIOLOGICAL_DEMAND.f26bp_uM
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `julia --project=. -e 'using Pkg; Pkg.activate("."); include("test/test_grid_io.jl")'`

Expected: FAIL. `ADP_uM` is not in the columns, and there is no `f26bp_uM` key.

- [ ] **Step 3: Gate on the PPP.jl commit, then start both regenerations in the background**

```bash
git -C /home/james/projects/PentosePhosphatePathway.jl rev-parse HEAD      # fcccae8a196ce9cfba23fae91ca77c4be4540807
git -C /home/james/projects/PentosePhosphatePathway.jl status --porcelain  # empty
SCRATCH=/tmp/claude-1000/-home-james-projects-PPP-Experiments/a3d4eca5-2f02-4011-8c30-18a62260a8ba/scratchpad/regen
mkdir -p $SCRATCH/f0 $SCRATCH/f5
cat > $SCRATCH/regen.jl <<'EOF'
using GlycolysisPentosePathwayPhaseMap, PentosePhosphatePathway
dir, f26 = ARGS[1], parse(Float64, ARGS[2])
regenerate_grid(; atp_out = joinpath(dir, "atp_demand_grid.csv"),
                  demand_out = joinpath(dir, "demand_grid.csv"), f26bp_uM = f26)
EOF
```

Run each with `run_in_background: true` and wait for the completion notifications. Do not poll.

```bash
julia --project=/home/james/projects/PPP_Experiments $SCRATCH/regen.jl $SCRATCH/f0 0.0 > $SCRATCH/f0.log 2>&1
julia --project=/home/james/projects/PPP_Experiments $SCRATCH/regen.jl $SCRATCH/f5 5.0 > $SCRATCH/f5.log 2>&1
```

Each takes minutes. They may run concurrently.

- [ ] **Step 4: Drift check — F26BP 0 must reproduce today's grids exactly**

```bash
cat > $SCRATCH/drift_check.jl <<'EOF'
using CSV, DataFrames
old_dir, new_dir = ARGS
ok = true
for f in ("atp_demand_grid.csv", "demand_grid.csv")
    old = CSV.read(joinpath(old_dir, f), DataFrame); new = CSV.read(joinpath(new_dir, f), DataFrame)
    nrow(old) == nrow(new) || (println("$f: row count $(nrow(old)) -> $(nrow(new))"); global ok = false)
    for c in names(old)
        c in names(new) || (println("$f: column $c missing"); global ok = false; continue)
        isequal(old[!, c], new[!, c]) && continue
        println("$f: column $c differs"); global ok = false
    end
    added = setdiff(names(new), names(old))
    println("$f: added columns ", added)
    sort(added) == sort(["ADP_uM", "AMP_uM", "NAD_uM", "NADH_uM"]) || (global ok = false)
end
println(ok ? "DRIFT CHECK PASSED" : "DRIFT CHECK FAILED")
EOF
julia --project=/home/james/projects/PPP_Experiments $SCRATCH/drift_check.jl data $SCRATCH/f0
```

Expected: `added columns ["ADP_uM", "AMP_uM", "NAD_uM", "NADH_uM"]` for both files, then `DRIFT CHECK PASSED`. Comparing by name also covers `i_nadph`/`i_r5p`, so row order is checked too. If it prints `FAILED`, **stop and report** the differing columns. Something other than this change moved the model.

- [ ] **Step 5: Summarise the F26BP 0 → 5 µM change**

```bash
cat > $SCRATCH/f26bp_summary.jl <<'EOF'
using CSV, DataFrames, GlycolysisPentosePathwayPhaseMap
const G = GlycolysisPentosePathwayPhaseMap
old = G.load_grid(ARGS[1]); new = G.load_grid(ARGS[2])
collapsed(df) = (df.retcode .== :Terminated) .& (df.ATP_uM .< 100)   # healthy ATP is ~7 mM
println("| ATP frac | collapsed F26BP 0 | collapsed F26BP 5 | non-converged 0 | non-converged 5 |")
println("|---|---|---|---|---|")
for a in sort(unique(old.atpase_frac))
    o = old[old.atpase_frac .== a, :]; n = new[new.atpase_frac .== a, :]
    println("| ", round(a, sigdigits = 3), " | ", count(collapsed(o)), " | ", count(collapsed(n)),
            " | ", count(!=(:Terminated), o.retcode), " | ", count(!=(:Terminated), n.retcode), " |")
end
key(df) = [(r.atpase_frac, r.i_nadph, r.i_r5p) for r in eachrow(df)]
@assert key(old) == key(new)
println("cells changing mode: ", count(old.mode .!= new.mode), " of ", nrow(old))
in_, ir = G.physiological_cell(new); a = G._PHYSIOLOGICAL_DEMAND.atpase_frac
pick(df) = only(eachrow(df[(df.atpase_frac .≈ a) .& (df.i_nadph .== in_) .& (df.i_r5p .== ir), :]))
po, pn = pick(old), pick(new)
println("physiological cell: cycle_index ", po.cycle_index, " -> ", pn.cycle_index,
        "; mode ", po.mode, " -> ", pn.mode, "; ADP ", pn.ADP_uM, " AMP ", pn.AMP_uM, " ATP ", pn.ATP_uM)
low = new[(new.atpase_frac .≈ a) .& (new.i_nadph .== 1), :]
println("lowest-NADPH-demand row at ATP $a: collapsed ", count(collapsed(low)), " of ", nrow(low))
for df in (new,), nm in (:NADP_uM, :R5P_uM)
    alive = df[(df.retcode .== :Terminated) .& (df.ATP_uM .>= 100), nm]
    println("min live $nm = ", minimum(alive), " (floor 1e-4)")
end
EOF
julia --project=. $SCRATCH/f26bp_summary.jl data/atp_demand_grid.csv $SCRATCH/f5/atp_demand_grid.csv | tee $SCRATCH/f26bp_summary.md
```

Expected: the physiological-cell ADP is about 29 µM and AMP about 0.15 µM. Keep the whole output; it goes into the commit body. If a live (non-collapsed) cell has NADP or R5P below 1e-4 µM, note it. The spec's floor rationale (§4) and Task 6's docstring must then say so.

- [ ] **Step 6: Install the ship-run grids and their `[grids]` table**

```bash
cp $SCRATCH/f5/atp_demand_grid.csv $SCRATCH/f5/demand_grid.csv data/
julia --project=. -e '
using GlycolysisPentosePathwayPhaseMap; const G = GlycolysisPentosePathwayPhaseMap
g = G.read_provenance(ARGS[1])["grids"]
G.update_provenance!(G.provenance_path(), "grids", g)' $SCRATCH/f5/provenance.toml
git diff --stat data/
```

Expected: both CSVs and `data/provenance.toml` change, and `[literature]` is unchanged in `git diff data/provenance.toml`.

- [ ] **Step 7: Run the full core suite**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`

Expected: PASS, including the new testset. `test_phasemap.jl`'s "most strongly reverse cell has mode reverse" and `test_grid_io.jl`'s "physiological cell alive" now run on the 5 µM grids. If either fails, **stop and report**; do not edit the assertion.

- [ ] **Step 8: Commit**

```bash
git branch --show-current
git add data/atp_demand_grid.csv data/demand_grid.csv data/provenance.toml test/test_grid_io.jl
git commit -F - <<EOF
data(grids): regenerate both grids at 5 uM F26BP with the cofactor pools

Plain English: the shipped phase-map grids are re-solved with F26BP at its physiological level
and now carry ADP, AMP, NAD+ and NADH. A first run at F26BP 0 reproduced the old grids exactly,
so every change below comes from F26BP alone.

PPP.jl main fcccae8, axes unchanged. Drift check at f26bp_uM = 0: all existing columns
bit-identical (isequal by name), only ADP_uM, AMP_uM, NAD_uM, NADH_uM added.

$(cat $SCRATCH/f26bp_summary.md)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 3: Panel tuples and grid column checks

**Repo:** PhaseMap, branch `feat/free-pool-panel-split`.

**Files:**
- Modify: `src/grid_io.jl:5-7` (the panel order) and `:110-118` (`load_grid`)
- Modify: `src/explorer.jl:32-37` (the start of `build_explorer`)
- Test: `test/test_grid_io.jl`, `test/test_metab_panel.jl:12`, `test/test_phasemap.jl`

**Interfaces:**
- Produces: `_FREE_EQ_TOTAL_PANEL`, `_FREE_POOL_PANEL`, `_PANEL_METAB_ORDER` (length 20).
- Produces: `_require_columns(df, cols, why) -> df`, which throws `ArgumentError`.
- Produces: `_GRID_CORE_COLUMNS::Tuple{Vararg{String}}`.

- [ ] **Step 1: Write the failing tests**

Append to `test/test_grid_io.jl`:

```julia
@testset "panel tuples partition the panel order" begin
    l, r = G._FREE_EQ_TOTAL_PANEL, G._FREE_POOL_PANEL
    @test r == (:ADP, :AMP, :NAD, :NADH, :NADP, :NADPH)
    @test length(l) == 14 && isempty(intersect(l, r))
    @test G._PANEL_METAB_ORDER == (l..., r...)
    for path in (G.default_atp_grid(), G.default_demand_grid())
        cols = names(CSV.read(path, DataFrame; limit = 1))
        @test all("$(s)_uM" in cols for s in G._PANEL_METAB_ORDER)
    end
end

@testset "grid column checks" begin
    df = CSV.read(G.default_atp_grid(), DataFrame)
    p = tempname() * ".csv"
    CSV.write(p, select(df, Not(:V_G6PD_uM_min)))
    err = try G.load_grid(p); nothing catch e; e end
    @test err isa ArgumentError && occursin("V_G6PD_uM_min", err.msg)

    p2 = tempname() * ".csv"
    CSV.write(p2, select(df, Not([:ADP_uM, :NADH_uM])))
    old = G.load_grid(p2)                                   # load_grid does not need pools
    err2 = try G.build_explorer(old); nothing catch e; e end
    @test err2 isa ArgumentError
    @test occursin("ADP_uM", err2.msg) && occursin("NADH_uM", err2.msg)
    @test occursin("regenerate", err2.msg)

    p3 = tempname() * ".csv"
    CSV.write(p3, select(df, Not(:atpase_frac)))
    err3 = try G.build_explorer(G.load_grid(p3)); nothing catch e; e end
    @test err3 isa ArgumentError && occursin("atpase_frac", err3.msg)
end
```

Append to `test/test_phasemap.jl`:

```julia
@testset "phasemap() renders a grid without pool columns" begin
    df = CSV.read(G.default_demand_grid(), DataFrame)
    p = tempname() * ".csv"
    CSV.write(p, select(df, Not([:ADP_uM, :AMP_uM, :NAD_uM, :NADH_uM])))
    @test phasemap(grid = p) isa Figure
end
```

Also change the first line of `test/test_phasemap.jl` to `using Test, CSV, DataFrames, CairoMakie, GlycolysisPentosePathwayPhaseMap`.

In `test/test_metab_panel.jl`, replace line 12, which fails once the panel holds 20 metabolites and the literature table doesn't yet have the cofactors:

```julia
    @test all(s -> !isempty(lit[s]), filter(!=(:PGLn), collect(G._FREE_EQ_TOTAL_PANEL)))
```

- [ ] **Step 2: Run them to verify they fail**

Run: `julia --project=. -e 'using Pkg; Pkg.activate("."); include("test/test_grid_io.jl")'`

Expected: FAIL with `UndefVarError: _FREE_EQ_TOTAL_PANEL`. `load_grid` does not throw.

- [ ] **Step 3: Implement the tuples**

In `src/grid_io.jl`, replace lines 5-7:

```julia
# Metabolite panel, split as in PPP_Experiments' basic panel. Left: species whose free and
# total pools match, compared with extraction metabolomics. Right: cofactors whose pool is
# mostly protein-bound, compared with free-pool data only. The grouping is PPP_Experiments'
# scripts/helper_files/metab_data.jl (FREE_EQ_TOTAL_METABS / FREE_POOL_METABS); its literature
# export (check_panel_grouping) refuses to run if these tuples disagree with it.
const _FREE_EQ_TOTAL_PANEL = (:Glucose, :G6P, :F6P, :F16BP, :DHAP, :GAP, :ATP,
                              :PGLn, :PGA, :Ru5P, :R5P, :X5P, :S7P, :E4P)
const _FREE_POOL_PANEL     = (:ADP, :AMP, :NAD, :NADH, :NADP, :NADPH)
# Display order = the two panels concatenated. Single source of truth for row_to_pools,
# load_metab_lit, the metabolite panel and the PPP_Experiments export.
const _PANEL_METAB_ORDER   = (_FREE_EQ_TOTAL_PANEL..., _FREE_POOL_PANEL...)
```

- [ ] **Step 4: Implement the column checks**

In `src/grid_io.jl`, replace `load_grid` and its docstring:

```julia
# Columns every renderer reads (phasemap and explore). Flux columns come from _PPP_FLUX_KEYS
# (flux_map.jl, included after this file; resolved at call time).
const _GRID_CORE_COLUMNS = ("retcode", "cycle_index", "mode", "i_nadph", "i_r5p",
                            "nadph_phi", "r5p_phi")

"Throw an `ArgumentError` naming every entry of `cols` missing from `df`, with reason `why`."
function _require_columns(df, cols, why)
    absent = [c for c in cols if !(c in names(df))]
    isempty(absent) || throw(ArgumentError("grid is missing column(s) $(join(absent, ", ")): $why"))
    return df
end

"""
    load_grid(path) -> DataFrame

Read a grid CSV, check it has the columns every renderer needs (`_GRID_CORE_COLUMNS` and a
`<key>_uM_min` flux column per `_PPP_FLUX_KEYS`), and normalize `retcode` from String to Symbol
(CSV round-trips it as "Terminated"; the renderers mask on the `:Terminated` Symbol). Pool
columns are checked by `build_explorer`, the only reader.
"""
function load_grid(path)
    df = CSV.read(path, DataFrame)
    _require_columns(df, [_GRID_CORE_COLUMNS..., ("$(k)_uM_min" for k in _PPP_FLUX_KEYS)...],
                     "not a phase-map grid CSV")
    df.retcode = Symbol.(df.retcode)
    return df
end
```

In `src/explorer.jl`, make these the first lines of `build_explorer`'s body:

```julia
    _require_columns(df, ["atpase_frac"],
                     "the explorer needs the 3-D ATP grid (atp_demand_grid.csv)")
    _require_columns(df, ["$(s)_uM" for s in _PANEL_METAB_ORDER],
                     "the grid predates the cofactor pool columns; regenerate it with regenerate_grid()")
```

- [ ] **Step 5: Run the full core suite**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`

Expected: PASS. The existing `draw_metab_panel!` still runs, now with 20 columns and empty boxes for the cofactors.

- [ ] **Step 6: Commit**

```bash
git branch --show-current
git add src/grid_io.jl src/explorer.jl test/test_grid_io.jl test/test_metab_panel.jl test/test_phasemap.jl
git commit -F - <<'EOF'
feat(grid-io): split the panel into free ~ total and free-pool tuples; check grid columns

Plain English: the explorer now knows which metabolites belong on the "free equals total"
side and which on the free-pool side. Opening an old grid that lacks the new cofactor columns
fails with a clear message to regenerate it; the static phase map still accepts old grids.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 4: PPP_Experiments — matched-pair export, `method` column, grouping guard

**Repo:** PPP_Experiments, branch `feat/phasemap-free-pool-export`. PhaseMap must be on `feat/free-pool-panel-split` with Task 3 committed: the environment `develop`s it by path, and the export reads its tuples.

**Files:**
- Modify: `scripts/helper_files/metab_data.jl:287-317` (`derive_matched_free`)
- Modify: `scripts/helper_files/test/test_free_pool.jl` (after line 97)
- Modify: `scripts/phasemap_export/literature_export.jl`, `scripts/phasemap_export/export_literature_pools.jl`
- Test: `scripts/phasemap_export/test_literature_export.jl`

**Interfaces:**
- Produces: `derive_matched_free` output gains columns `Ratio::String` (the ratio's `Quantity`, e.g. `"NADPH/NADP"`) and `Partner::String` (the measured side, e.g. `"NADPH"`). Existing columns are unchanged.
- Produces: `check_panel_grouping(free_eq_total_panel, free_pool_panel) -> nothing`, which throws `ArgumentError`.
- Produces: `literature_pools_table(metab_df, free_df, free_derived, panel) -> DataFrame` with columns `metabolite, value_uM, reference, cell_type, derived, method`.
- Produces: `derived_reference(r) -> String`.

- [ ] **Step 1: Write the failing `derive_matched_free` test**

In `scripts/helper_files/test/test_free_pool.jl`, add after the "partner is the numerator" testset:

```julia
@testset "derive_matched_free: records the ratio and the partner" begin
    d = _derive([("Tao", "HeLa", "cytosol", "NADPH", 3.1, "uM", "sensor"),
                 ("Sallin", "HeLa", "cytosol", "NADPH/NADP", 75.0, "ratio", "sensor")])
    @test d.Ratio[1] == "NADPH/NADP" && d.Partner[1] == "NADPH"
    d2 = _derive([("X", "HeLa", "cytosol", "NADP", 0.04, "uM", "sensor"),
                  ("Sallin", "HeLa", "cytosol", "NADPH/NADP", 75.0, "ratio", "sensor")])
    @test d2.Partner[1] == "NADP"
end
```

- [ ] **Step 2: Rewrite the export tests (they fail against the current export)**

Replace everything in `scripts/phasemap_export/test_literature_export.jl` from `# free_df as prepare_free_pool leaves it` to the end of the file with the following. Keep the header, the `using` line (with `Random` removed), the two `include`s and `metab_df()`.

```julia
# free_df as prepare_free_pool leaves it: concentrations in M.
function free_df(; extra = DataFrame())
    df = DataFrame("Reference" => ["Tao", "Moon", "Mito", "Sallin 4", "Sallin NAD"],
                   "Cell Type" => ["HeLa", "HeLa", "HeLa", "HeLa", "HeLa"],
                   "Compartment" => ["cytosol", "cytosol", "mitochondria", "cytosol", "cytosol"],
                   "Quantity" => ["NADPH", "NADPH", "NADPH", "NADPH/NADP", "NAD"],
                   "Value" => [3.1e-6, 2.6e-6, 9.9e-6, 75.0, 70e-6],
                   "Unit" => ["M", "M", "M", "ratio", "M"],
                   "Method" => ["sensor", "sensor", "sensor", "sensor", "sensor"])
    return vcat(df, extra; cols = :union)
end

const PANEL = (:G6P, :PGLn, :NADPH, :NADP, :AMP, :NAD)

@testset "literature_pools_table" begin
    fdf = free_df()
    t = literature_pools_table(metab_df(), fdf, derive_matched_free(fdf, metab_df()), PANEL)
    @test names(t) == ["metabolite", "value_uM", "reference", "cell_type", "derived", "method"]

    g = t[t.metabolite .== "G6P", :]
    @test g.value_uM ≈ [100.0, 50.0]                          # M -> µM; missing dropped
    @test g.cell_type == ["HeLa", "HeLa_Oligomycin"]          # treated rows kept
    @test !any(g.derived) && all(==("metabolomics"), g.method)
    @test !any(t.metabolite .== "PGLn")                       # no data -> no rows
    @test !any(t.metabolite .== "AMP")                        # free-pool, no data -> no rows, no error

    n = t[t.metabolite .== "NADPH", :]
    @test n.value_uM ≈ [3.1, 2.6]                             # cytosol only (mito dropped)
    @test n.reference == ["Tao", "Moon"] && all(==("sensor"), n.method) && !any(n.derived)

    d = t[t.metabolite .== "NADP", :]
    @test nrow(d) == 1 && d.derived[1] && d.method[1] == "sensor" && d.cell_type[1] == "HeLa"
    @test d.value_uM[1] ≈ median([3.1, 2.6]) / 75.0
    @test d.reference[1] == "derived: NADPH/NADP ratio from Sallin 4; median NADPH of n=2 same-line rows"

    @test all(==("sensor"), t.method[t.metabolite .== "NAD"])
    path = tempname() * ".csv"
    CSV.write(path, t)
    @test isascii(read(path, String))
end

@testset "literature_pools_table guards" begin
    fdf = free_df()
    fd = derive_matched_free(fdf, metab_df())
    @test_throws ArgumentError literature_pools_table(metab_df(), fdf, fd, (:Foo,))
    bad = metab_df(); bad.G6P[1] = 0.0
    @test_throws ArgumentError literature_pools_table(bad, fdf, fd, (:G6P,))
    nonascii = metab_df(); nonascii.Reference[1] = "Müller"
    @test_throws ArgumentError literature_pools_table(nonascii, fdf, fd, (:G6P,))
end

@testset "check_panel_grouping" begin
    @test check_panel_grouping((:G6P, :ATP), FREE_POOL_METABS) === nothing
    @test_throws ArgumentError check_panel_grouping((:G6P,), reverse(FREE_POOL_METABS))
    @test_throws ArgumentError check_panel_grouping((:G6P,), FREE_POOL_METABS[1:5])
    @test_throws ArgumentError check_panel_grouping((:G6P, :NADH), FREE_POOL_METABS)
    @test_throws ArgumentError check_panel_grouping((:Foo,), FREE_POOL_METABS)
end
```

- [ ] **Step 3: Run both tests to verify they fail**

```bash
cd /home/james/projects/PPP_Experiments
julia --project scripts/helper_files/test/test_free_pool.jl
julia --project scripts/phasemap_export/test_literature_export.jl
```

Expected: FAIL. The first fails with no `Ratio`/`Partner` column. The second fails with a `MethodError` (four-argument `literature_pools_table`) and an undefined `check_panel_grouping`.

- [ ] **Step 4: Add `Ratio` and `Partner` to `derive_matched_free`**

In `scripts/helper_files/metab_data.jl`, change the docstring's column list to "Columns: Quantity, Value (M), Cell Type, Method (the ratio's), Reference (the ratio's), n_partner_rows, Ratio (the ratio's Quantity), Partner (the measured side)." In the body:

```julia
    out = DataFrame("Quantity" => String[], "Value" => Float64[], "Cell Type" => String[],
                    "Method" => String[], "Reference" => String[], "n_partner_rows" => Int[],
                    "Ratio" => String[], "Partner" => String[])
```

and replace the two lines that compute and push the derived row:

```julia
        derived, partner, value, n = isempty(pd) ? (den, num, median(pn) / r.Value, length(pn)) :
                                                   (num, den, median(pd) * r.Value, length(pd))
        push!(out, (string(derived), value, cell, r.Method, r.Reference, n,
                    String(r.Quantity), string(partner)))
```

- [ ] **Step 5: Rewrite `literature_export.jl`**

Replace the file's body below the banner. Keep the banner, and update its header comment bullets to match.

```julia
# Builds the long-format table GlycolysisPentosePathwayPhaseMap.jl ships as
# data/literature_pools.csv: one row per literature value, in µM on the model's pool basis.
# Inputs are the frames import_data.jl leaves behind (values in M, already converted by row
# Basis and species class, GAP/DHAP free-aldehyde corrected), so the basis conversion lives
# only in metab_data.jl. Caller must have included metab_data.jl.
#
# - FREE_EQ_TOTAL_METABS on the panel: every metabolomics row, treated conditions included.
# - FREE_POOL_METABS on the panel: direct cytosolic free-pool rows, then the matched-pair
#   derived rows (derive_matched_free); never the total-pool metabolomics columns.
# - `method` tells the PhaseMap viewer which marker to draw: metabolomics, sensor, equilibrium.

"Cytosolic free-pool concentration rows (M) of metabolite `m`."
_free_rows(free_df, m::Symbol) = free_df[(free_df.Compartment .== "cytosol") .&
                                         (free_df.Unit .== "M") .&
                                         (free_df.Quantity .== string(m)), :]

"Reference text of a matched-pair derived row: both inputs, the ratio and its partner."
derived_reference(r) =
    "derived: $(r.Ratio) ratio from $(r.Reference); median $(r.Partner) of n=$(r.n_partner_rows) same-line rows"

"""
    check_panel_grouping(free_eq_total_panel, free_pool_panel) -> nothing

Throw an `ArgumentError` unless PhaseMap's panel grouping agrees with metab_data.jl:
`free_pool_panel` equals `FREE_POOL_METABS` (same order), and every `free_eq_total_panel`
entry is in `FREE_EQ_TOTAL_METABS`.
"""
function check_panel_grouping(free_eq_total_panel, free_pool_panel)
    Tuple(free_pool_panel) == FREE_POOL_METABS || throw(ArgumentError(
        "PhaseMap _FREE_POOL_PANEL $(Tuple(free_pool_panel)) != FREE_POOL_METABS " *
        "$(FREE_POOL_METABS) (metab_data.jl); make them match before exporting"))
    bad = [m for m in free_eq_total_panel if !(m in FREE_EQ_TOTAL_METABS)]
    isempty(bad) || throw(ArgumentError(
        "PhaseMap _FREE_EQ_TOTAL_PANEL entries not in FREE_EQ_TOTAL_METABS (metab_data.jl): $bad"))
    return nothing
end

"""
    literature_pools_table(metab_df, free_df, free_derived, panel) -> DataFrame

Rows `metabolite, value_uM, reference, cell_type, derived, method` for every metabolite in
`panel`. `free_derived` is `derive_matched_free(free_df, metab_df)` (import_data.jl's
`free_pool_derived`). A free-pool metabolite with no rows contributes none. Throws
`ArgumentError` when a panel metabolite is neither a free-pool metabolite nor a metabolomics
column, when any value is not finite and > 0, or when any text field is not ASCII.
"""
function literature_pools_table(metab_df, free_df, free_derived, panel)
    out = DataFrame(metabolite = String[], value_uM = Float64[], reference = String[],
                    cell_type = String[], derived = Bool[], method = String[])
    mcols = metab_columns(metab_df)
    for m in panel
        if m in FREE_POOL_METABS
            for r in eachrow(_free_rows(free_df, m))
                push!(out, (string(m), r.Value * 1e6, String(r.Reference),
                            String(r["Cell Type"]), false, String(r.Method)))
            end
            for r in eachrow(free_derived[free_derived.Quantity .== string(m), :])
                push!(out, (string(m), r.Value * 1e6, derived_reference(r),
                            String(r["Cell Type"]), true, String(r.Method)))
            end
        elseif string(m) in mcols
            for r in eachrow(metab_df)
                v = r[string(m)]
                ismissing(v) && continue
                push!(out, (string(m), v * 1e6, String(r.Reference), String(r["Cell Type"]),
                            false, "metabolomics"))
            end
        else
            throw(ArgumentError("panel metabolite :$m has no literature source (no " *
                                "metabolomics column and not a free-pool metabolite)"))
        end
    end
    bad = findall(v -> !(isfinite(v) && v > 0), out.value_uM)
    isempty(bad) || throw(ArgumentError(
        "non-finite or non-positive literature values: $(collect(eachrow(out[bad, :])))"))
    nonascii = [s for s in vcat(out.reference, out.cell_type) if !isascii(s)]
    isempty(nonascii) || throw(ArgumentError("non-ASCII text in literature table: $(unique(nonascii))"))
    return out
end
```

- [ ] **Step 6: Update `export_literature_pools.jl`**

- Change `using PentosePhosphatePathway, CSV, DataFrames, Random` to `using PentosePhosphatePathway, CSV, DataFrames`.
- Delete the `NADP_DRAWS` and `NADP_SEED` constants.
- Replace the `tbl = literature_pools_table(...)` line with:

```julia
check_panel_grouping(GPPM._FREE_EQ_TOTAL_PANEL, GPPM._FREE_POOL_PANEL)
tbl = literature_pools_table(metab_conc_data, free_pool_data, free_pool_derived,
                             GPPM._PANEL_METAB_ORDER)
```

- In the provenance `Dict`, delete the `"nadp_seed"` and `"nadp_draws"` entries.
- Update the header comment: the import line now reads `# metab_conc_data, free_pool_data, free_pool_derived`.

- [ ] **Step 7: Run the tests to verify they pass**

```bash
julia --project scripts/helper_files/test/test_free_pool.jl
julia --project scripts/phasemap_export/test_literature_export.jl
julia --project scripts/helper_files/test/runtests.jl
```

Expected: all PASS. `runtests.jl` covers `test_plot_metab_split.jl`, which also consumes `derive_matched_free`.

- [ ] **Step 8: Dry-run the export into scratch and check the result**

```bash
SCRATCH=/tmp/claude-1000/-home-james-projects-PPP-Experiments/a3d4eca5-2f02-4011-8c30-18a62260a8ba/scratchpad/regen
mkdir -p $SCRATCH/lit && cp ../GlycolysisPentosePathwayPhaseMap.jl/data/provenance.toml $SCRATCH/lit/
julia --project scripts/phasemap_export/export_literature_pools.jl $SCRATCH/lit
```

Expected: `wrote N rows`. The per-metabolite counts must include ADP 6, NAD 6, NADH 2, NADP 1 and NADPH 2, and no AMP row. This is a dry run only; the shipped file is written in Task 5.

- [ ] **Step 9: Commit, push and open the PR**

```bash
git branch --show-current   # feat/phasemap-free-pool-export
git add scripts/helper_files/metab_data.jl scripts/helper_files/test/test_free_pool.jl scripts/phasemap_export/literature_export.jl scripts/phasemap_export/export_literature_pools.jl scripts/phasemap_export/test_literature_export.jl
git commit -F - <<'EOF'
feat(phasemap-export): export free-pool cofactors as matched-pair points with a method column

Plain English: the literature table shipped to the phase-map viewer now covers all six
protein-bound cofactors from free-pool measurements only, and says for every value whether it
came from metabolomics, a sensor or an equilibrium estimate, so the viewer can mark them like
the basic panel does. NADP becomes the single matched HeLa point instead of 1000 random draws.
The export refuses to run if the viewer's grouping differs from this repo's.

derive_matched_free gains Ratio and Partner columns (additive). literature_pools_table takes
free_derived; NADP draw/seed machinery and its provenance keys removed; check_panel_grouping
guards GPPM._FREE_EQ_TOTAL_PANEL / _FREE_POOL_PANEL against metab_data.jl.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
git push -u origin feat/phasemap-free-pool-export
gh pr create --title "feat(phasemap-export): free-pool cofactor rows and method column" --body "<plain-English summary + the commit body above>

🤖 Generated with [Claude Code](https://claude.com/claude-code)"
```

- [ ] **Step 10: Checkpoint — the user merges the PR**

Ask the user to merge the PR (squash), or to confirm that you may. Do not continue to Task 5 until it is merged.

---

### Task 5: Ship the literature table from PPP_Experiments `main`

**Repo:** run in PPP_Experiments `main`; commit to PhaseMap `feat/free-pool-panel-split`.

**Files:**
- Replace: PhaseMap `data/literature_pools.csv`
- Modify: PhaseMap `data/provenance.toml` (`[literature]` only)
- Modify: PhaseMap `src/metab_panel.jl` (the `_METAB_YLIMS` constant only)
- Test: PhaseMap `test/test_metab_panel.jl:5-17`

**Interfaces:**
- Produces: a shipped table with a `method` column and free-pool rows; the constant `_METAB_YLIMS = (1e-4, 1e5)`.

- [ ] **Step 1: Write the failing shipped-table tests**

In PhaseMap `test/test_metab_panel.jl`, replace the body of the testset `"load_metab_lit (shipped literature_pools.csv)"` from `raw = CSV.read(...)` through `@test !any(raw.derived[raw.metabolite .!= "NADP"])` with:

```julia
    raw = CSV.read(path, DataFrame)
    @test names(raw) == ["metabolite", "value_uM", "reference", "cell_type", "derived", "method"]
    @test Set(raw.method) ⊆ Set(["metabolomics", "sensor", "equilibrium"])
    left  = in(string.(G._FREE_EQ_TOTAL_PANEL)).(raw.metabolite)
    right = in(string.(G._FREE_POOL_PANEL)).(raw.metabolite)
    @test all(raw.method[left] .== "metabolomics") && !any(raw.derived[left])
    @test !any(raw.method[right] .== "metabolomics")
    @test all(raw.derived[raw.metabolite .== "NADP"]) && any(raw.metabolite .== "NADP")
    @test !any(raw.metabolite .== "AMP")
    @test all(s -> any(raw.metabolite .== string(s)), (:ADP, :NAD, :NADH, :NADPH))
    lo, hi = G._METAB_YLIMS                                   # no literature point is off-scale
    @test all(v -> lo <= v <= hi, raw.value_uM)
```

Also replace the earlier line `@test all(s -> !isempty(lit[s]), filter(!=(:PGLn), collect(G._FREE_EQ_TOTAL_PANEL)))` with:

```julia
    @test all(s -> !isempty(lit[s]), filter(∉((:PGLn, :AMP)), collect(G._PANEL_METAB_ORDER)))
```

- [ ] **Step 2: Run it to verify it fails**

Run: `julia --project=. -e 'using Pkg; Pkg.activate("."); include("test/test_metab_panel.jl")'` (in PhaseMap).

Expected: FAIL. There is no `method` column and `_METAB_YLIMS` is undefined.

- [ ] **Step 3: Add `_METAB_YLIMS`**

In PhaseMap `src/metab_panel.jl`, after `_METAB_YLABEL`:

```julia
"Metabolite-panel y-limits, µM. 1e-4 keeps every non-collapsed grid cell's pools on scale."
const _METAB_YLIMS = (1e-4, 1e5)
```

and change `ylims!(ax, 1e-3, 1e5)` to `ylims!(ax, _METAB_YLIMS...)`.

- [ ] **Step 4: Export from the merged squash commit**

```bash
cd /home/james/projects/PPP_Experiments
git switch main && git pull --ff-only
git log --oneline -1          # the PR's squash commit
git status --porcelain        # empty, so [literature] records dirty = false
julia --project scripts/phasemap_export/export_literature_pools.jl
```

Expected: `wrote N rows` into `../GlycolysisPentosePathwayPhaseMap.jl/data/`. Then, in PhaseMap:

```bash
cd /home/james/projects/GlycolysisPentosePathwayPhaseMap.jl
grep -A12 '^\[literature\]' data/provenance.toml   # commit = the squash hash, dirty = false, branch = main, no nadp_* keys
```

- [ ] **Step 5: Run the full core suite**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git branch --show-current   # feat/free-pool-panel-split
git add data/literature_pools.csv data/provenance.toml src/metab_panel.jl test/test_metab_panel.jl
git commit -F - <<'EOF'
data(literature): ship free-pool cofactor rows with a method column

Plain English: the viewer's literature table now holds free-pool measurements for ADP, NAD+,
NADH, NADP+ and NADPH (AMP has none), each tagged as metabolomics, sensor or equilibrium.
Exported from PPP_Experiments main at <squash hash> (PR #<n>).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

Replace `<squash hash>` and `<n>` with the real values before committing.

---

### Task 6: Two-axis metabolite panel

**Repo:** PhaseMap, branch `feat/free-pool-panel-split`.

**Files:**
- Modify: `src/metab_panel.jl` (rewrite below the constants)
- Modify: `src/explorer.jl:74-83`
- Test: `test/test_metab_panel.jl`, `test/test_explorer.jl:30-33`

**Interfaces:**
- Consumes: `_FREE_EQ_TOTAL_PANEL`, `_FREE_POOL_PANEL`, `_PANEL_METAB_ORDER` (Task 3); `_METAB_YLIMS` (Task 5).
- Produces:
  - `load_metab_lit(path; order) -> Dict{Symbol, Vector{_LitPoint}}`, where `_LitPoint = @NamedTuple{value::Float64, derived::Bool, method::Symbol}`;
  - `_free_pool_point_style(derived::Bool, method::Symbol) -> (; label, marker, color, strokecolor, strokewidth)`;
  - `_offscale_marks(pools; lims) -> (; on, below, above, labels, label_pos)`;
  - `draw_metab_split!(gl, lit, metab_obs) -> (; ax_left, ax_right)`.
- Removes: `draw_metab_panel!`.

- [ ] **Step 1: Write the failing tests**

In `test/test_metab_panel.jl`:

In the shipped testset, replace the two lines that assume `Vector{Float64}` values (today lines 10-11, `@test all(v -> v isa Vector{Float64}, values(lit))` and `@test all(v -> all(x -> isfinite(x) && x > 0, v), values(lit))`) with:

```julia
    @test all(v -> v isa Vector{G._LitPoint}, values(lit))
    @test all(v -> all(p -> isfinite(p.value) && p.value > 0, v), values(lit))
```

Then add right after the `raw` assertions in the same testset:

```julia
    @test all(p -> p.method === :metabolomics, lit[:G6P])
    @test only(lit[:NADP]).derived
```

Replace the testset `"load_metab_lit long format (synthetic)"` with:

```julia
@testset "load_metab_lit long format (synthetic)" begin
    path = tempname() * ".csv"
    CSV.write(path, DataFrame(metabolite = ["G6P", "G6P", "NADH", "Foo"], value_uM = [1.0, 2.0, 0.3, 3.0],
                              reference = ["a", "b", "c", "d"], cell_type = ["x", "y", "z", "w"],
                              derived = [false, false, true, false],
                              method = ["metabolomics", "metabolomics", "equilibrium", "sensor"]))
    lit = G.load_metab_lit(path)
    @test [p.value for p in lit[:G6P]] == [1.0, 2.0]   # no conversion applied by the loader
    @test only(lit[:NADH]) == (value = 0.3, derived = true, method = :equilibrium)
    @test isempty(lit[:PGLn]) && isempty(lit[:AMP])
    @test !haskey(lit, :Foo)

    wrong(m, meth) = (p = tempname() * ".csv";
        CSV.write(p, DataFrame(metabolite = [m], value_uM = [1.0], reference = ["a"],
                               cell_type = ["x"], derived = [false], method = [meth])); p)
    @test_throws ArgumentError G.load_metab_lit(wrong("G6P", "sensor"))       # left, not metabolomics
    @test_throws ArgumentError G.load_metab_lit(wrong("NADH", "metabolomics")) # total pool on the right
    @test_throws ArgumentError G.load_metab_lit(wrong("G6P", "guess"))        # unknown method
end

@testset "_free_pool_point_style (basic-panel marker key)" begin
    s = G._free_pool_point_style(false, :sensor)
    @test s.marker == :circle && s.color == :black
    @test G._free_pool_point_style(true, :sensor).color == :transparent
    @test G._free_pool_point_style(false, :equilibrium).marker == :diamond
end

@testset "_offscale_marks" begin
    m = G._offscale_marks([10.0, 1e-6, 0.0, -1e-18, NaN, 1e7]; lims = (1e-4, 1e5))
    @test m.on == [Point2f(1, 10.0)]
    @test [p[1] for p in m.below] == [2, 3, 4]               # tiny, zero and negative: floor
    @test [p[1] for p in m.above] == [6]                     # NaN draws nothing
    @test length(m.labels) == length(m.label_pos) == 4
    @test m.labels[2] == "0.0"
    @test all(p -> p[2] ≈ 1e-4 * 1.5, m.below) && all(p -> p[2] ≈ 1e5 / 1.5, m.above)
end
```

Replace the testset `"draw_metab_panel! offscreen"` with:

```julia
@testset "draw_metab_split! offscreen" begin
    lit = Dict(s => G._LitPoint[] for s in G._PANEL_METAB_ORDER)
    lit[:G6P] = [(value = v, derived = false, method = :metabolomics) for v in (30.0, 40.0, 55.0)]
    lit[:NADH] = [(value = 0.13, derived = false, method = :sensor),
                  (value = 0.34, derived = true, method = :equilibrium)]
    n = length(G._PANEL_METAB_ORDER)
    pools = fill(10.0, n); pools[2] = 1e-7; pools[16] = 0.0; pools[17] = NaN   # G6P off, AMP 0, NAD NaN
    metab_obs = Observable(pools)
    fig = Figure(size = (1500, 400)); gl = GridLayout(fig[1, 1])
    axes = G.draw_metab_split!(gl, lit, metab_obs)
    @test axes.ax_left.yscale[] === log10 && axes.ax_right.yscale[] === log10
    @test axes.ax_left.ylabel[] == G._METAB_YLABEL == "[metabolite], µM (cytosolic water)"
    @test occursin("no data", axes.ax_right.xticks[][2][2])                   # AMP column
    out = joinpath(tempdir(), "metab_split_smoke.png"); @test_nowarn save(out, fig); @test isfile(out)

    nplots = (length(axes.ax_left.scene.plots), length(axes.ax_right.scene.plots))
    for k in (3, 0, 5)                                        # off-scale count changes per update
        p = fill(10.0, n); p[1:k] .= 1e-9; p[end-k+1:end] .= 1e9
        @test_nowarn (metab_obs[] = p)
        @test_nowarn save(out, fig)
        @test (length(axes.ax_left.scene.plots), length(axes.ax_right.scene.plots)) == nplots
    end
    @test axes.ax_right.finallimits[].origin[1] ≈ 0.5
    @test axes.ax_right.finallimits[].widths[1] ≈ length(G._FREE_POOL_PANEL) + 0.3
    @test axes.ax_left.finallimits[].widths[1] ≈ length(G._FREE_EQ_TOTAL_PANEL) + 0.3
end
```

In `test/test_explorer.jl`, replace the `ax_metab = only(...)` block (lines 30-33) with:

```julia
    metab_axes = filter(x -> x isa Axis && x.yscale[] === log10 &&
                             (x.ylabel[] == G._METAB_YLABEL || occursin("free pool", x.title[])),
                        fig.content)
    @test length(metab_axes) == 2
    for ax in metab_axes, itr in (:rectanglezoom, :dragpan, :scrollzoom)
        @test !haskey(ax.interactions, itr)
    end
```

If `fig.content` does not list the nested axes, collect them with `filter(x -> x isa Axis, Makie.contents(fig.layout; exact = false))` instead. The review found that nested axes do appear in `fig.content`.

- [ ] **Step 2: Run them to verify they fail**

Run: `julia --project=. -e 'using Pkg; Pkg.activate("."); include("test/test_metab_panel.jl")'`

Expected: FAIL with `UndefVarError` for `_LitPoint`, `_free_pool_point_style`, `_offscale_marks` and `draw_metab_split!`.

- [ ] **Step 3: Rewrite `src/metab_panel.jl`**

Replace the whole file with:

```julia
# metab_panel.jl — metabolite panel: free ≈ total axis + free-pool axis, dynamic model marks

# The literature table is exported by PPP_Experiments (scripts/phasemap_export/) already on the
# model's pool basis: µM of cytosolic water, GAP/DHAP as the free aldehyde, and free-pool data
# only for the _FREE_POOL_PANEL cofactors. This package does no basis conversion. The split
# mirrors PPP_Experiments' basic panel (scripts/helper_files/plot_metab_split.jl).

"Metabolite-panel y-axis label; states the basis both the data and the model marks are on."
const _METAB_YLABEL = "[metabolite], µM (cytosolic water)"

"Metabolite-panel y-limits, µM. 1e-4 keeps every non-collapsed grid cell's pools on scale."
const _METAB_YLIMS = (1e-4, 1e5)

"Allowed `method` values in literature_pools.csv."
const _METAB_METHODS = (:metabolomics, :sensor, :equilibrium)

"Extra x-room right of the last column, for an off-scale value label."
const _XPAD_RIGHT = 0.3

"One literature value: µM, whether it was derived (matched pair), and how it was measured."
const _LitPoint = @NamedTuple{value::Float64, derived::Bool, method::Symbol}

"""
    load_metab_lit(csv_path; order = _PANEL_METAB_ORDER) -> Dict{Symbol, Vector{_LitPoint}}

Read the long-format literature table (`metabolite, value_uM, reference, cell_type, derived,
method`) and return, per panel metabolite, its values as `_LitPoint`s. A panel metabolite with
no rows maps to an empty vector; rows outside `order` are ignored. Throws `ArgumentError` on a
`method` outside `_METAB_METHODS`, on a non-metabolomics row for a `_FREE_EQ_TOTAL_PANEL`
metabolite, and on a metabolomics (total-pool) row for a `_FREE_POOL_PANEL` metabolite.
"""
function load_metab_lit(csv_path; order = _PANEL_METAB_ORDER)
    df = CSV.read(csv_path, DataFrame; types = Dict(:metabolite => String, :value_uM => Float64,
                                                    :derived => Bool, :method => String))
    lit = Dict{Symbol, Vector{_LitPoint}}(sym => _LitPoint[] for sym in order)
    for r in eachrow(df)
        sym, m = Symbol(r.metabolite), Symbol(r.method)
        m in _METAB_METHODS || throw(ArgumentError(
            "literature row for $sym has method :$m; expected one of $(_METAB_METHODS)"))
        haskey(lit, sym) || continue
        if sym in _FREE_POOL_PANEL
            m === :metabolomics && throw(ArgumentError(
                "free-pool metabolite $sym has a metabolomics (total-pool) row"))
        elseif m !== :metabolomics
            throw(ArgumentError("free ≈ total metabolite $sym has a non-metabolomics row (:$m)"))
        end
        push!(lit[sym], (value = r.value_uM, derived = r.derived, method = m))
    end
    return lit
end

"""
    _free_pool_point_style(derived, method) -> NamedTuple

Marker for one free-pool data point: filled = measured, open = inferred (by us or the original
authors); circle = sensor-based, diamond = equilibrium-based. Copied from PPP_Experiments
scripts/helper_files/plot_metab_split.jl (`free_pool_point_style`) so both panels share one key.
"""
function _free_pool_point_style(derived::Bool, method::Symbol)
    method === :equilibrium && return (label = "equilibrium-based", marker = :diamond,
        color = :transparent, strokecolor = :black, strokewidth = 1.5)
    derived && return (label = "derived free (sensor ratio)", marker = :circle,
        color = :transparent, strokecolor = :black, strokewidth = 1.5)
    return (label = "measured free (sensor)", marker = :circle, color = :black,
        strokecolor = :black, strokewidth = 0.0)
end

"""
    _offscale_marks(pools; lims = _METAB_YLIMS) -> (; on, below, above, labels, label_pos)

Split one axis's model pools (µM, column i = position i) into in-range diamond positions and
floor/ceiling arrow positions with value labels. Zero and negative solver values count as below
the floor; non-finite values (solver-error rows) draw nothing.
"""
function _offscale_marks(pools; lims = _METAB_YLIMS)
    lo, hi = lims
    on, below, above, label_pos = Point2f[], Point2f[], Point2f[], Point2f[]
    labels = String[]
    for (i, v) in enumerate(pools)
        isfinite(v) || continue
        if lo <= v <= hi
            push!(on, Point2f(i, v)); continue
        end
        y = v < lo ? lo * 1.5 : hi / 1.5
        push!(v < lo ? below : above, Point2f(i, y))
        push!(labels, string(round(v; sigdigits = 2))); push!(label_pos, Point2f(i + 0.12, y))
    end
    return (; on, below, above, labels, label_pos)
end

# One metabolite axis at gl[1, col]: categorical x pinned to the columns (sparse reactive marks
# must never auto-limit it), log y with the limits set BEFORE the log switch.
function _metab_axis(gl, col, title, labels)
    ax = Axis(gl[1, col]; title, xticks = (collect(eachindex(labels)), labels),
              xticklabelrotation = pi / 3)
    ylims!(ax, _METAB_YLIMS...)
    ax.yscale = log10
    xlims!(ax, 0.5, length(labels) + 0.5 + _XPAD_RIGHT)
    return ax
end

# Model marks lifted off `pools_obs`: a fixed set of four plots whose data change on reselection.
function _draw_model_marks!(ax, pools_obs)
    marks = lift(_offscale_marks, pools_obs)
    scatter!(ax, lift(m -> m.on, marks); marker = :diamond, markersize = 21, color = :purple,
             strokecolor = :black, strokewidth = 0.5)
    scatter!(ax, lift(m -> m.below, marks); marker = :dtriangle, markersize = 14, color = :purple)
    scatter!(ax, lift(m -> m.above, marks); marker = :utriangle, markersize = 14, color = :purple)
    text!(ax, lift(m -> m.label_pos, marks); text = lift(m -> m.labels, marks), fontsize = 10,
          color = :purple, align = (:left, :center))
end

# Static metabolomics boxes, one boxplot! call grouped by integer x-position.
function _draw_lit_boxes!(ax, lit, metabs)
    box_x = Float64[]; box_y = Float64[]
    for (i, sym) in enumerate(metabs), p in lit[sym]
        push!(box_x, i); push!(box_y, p.value)
    end
    isempty(box_x) || boxplot!(ax, box_x, box_y; width = 0.6, color = (:gray70, 0.6),
                               strokecolor = :gray30, strokewidth = 1, whiskerwidth = 0.6,
                               markersize = 4)
end

# Static free-pool points, spread over ±0.2 of the column when there are several.
function _draw_free_pool_points!(ax, lit, metabs)
    for (x, sym) in enumerate(metabs)
        pts = lit[sym]
        isempty(pts) && continue
        xs = length(pts) == 1 ? [Float64(x)] : collect(range(x - 0.2, x + 0.2, length = length(pts)))
        for (xi, p) in zip(xs, pts)
            s = _free_pool_point_style(p.derived, p.method)
            scatter!(ax, [xi], [p.value]; marker = s.marker, color = s.color,
                     strokecolor = s.strokecolor, strokewidth = s.strokewidth, markersize = 11)
        end
    end
end

function _metab_legend!(pos)
    styles = [_free_pool_point_style(false, :sensor), _free_pool_point_style(true, :sensor),
              _free_pool_point_style(false, :equilibrium)]
    elems = [PolyElement(color = (:gray70, 0.6), strokecolor = :gray30),
             MarkerElement(marker = :diamond, color = :purple, strokecolor = :black,
                           strokewidth = 0.5, markersize = 16),
             MarkerElement(marker = :dtriangle, color = :purple, markersize = 12),
             (MarkerElement(marker = s.marker, color = s.color, strokecolor = s.strokecolor,
                            strokewidth = s.strokewidth, markersize = 11) for s in styles)...]
    labels = ["literature (metabolomics)", "selected cell", "selected cell, off scale",
              (s.label for s in styles)...]
    return Legend(pos, elems, labels; framevisible = false, labelsize = 11, tellheight = false)
end

"""
    draw_metab_split!(gl, lit, metab_obs) -> (; ax_left, ax_right)

Render the metabolite panel into the `GridLayout` `gl`: left axis `_FREE_EQ_TOTAL_PANEL`
(static metabolomics boxes), right axis `_FREE_POOL_PANEL` (static free-pool points; a column
with none reads "no data"), legend in `gl[1, 3]` outside both data areas. `metab_obs` is a
length-`_PANEL_METAB_ORDER` Observable of the selected cell's pools (µM); its purple diamonds,
and labelled floor/ceiling triangles for off-scale values, update on reselection without
adding or removing plots. Linked log y over `_METAB_YLIMS`; widths follow metabolite counts.
"""
function draw_metab_split!(gl, lit, metab_obs)
    left, right = collect(_FREE_EQ_TOTAL_PANEL), collect(_FREE_POOL_PANEL)
    nl = length(left)
    right_labels = [isempty(lit[s]) ? "$s\nno data" : string(s) for s in right]
    ax_left = _metab_axis(gl, 1, "free ≈ total (extraction metabolomics)", string.(left))
    ax_left.ylabel = _METAB_YLABEL
    ax_right = _metab_axis(gl, 2, "free pool (free-pool data only)", right_labels)
    linkyaxes!(ax_left, ax_right)
    hideydecorations!(ax_right; grid = false, ticks = false)
    colsize!(gl, 1, Auto(length(left))); colsize!(gl, 2, Auto(length(right))); colgap!(gl, 12)
    _draw_lit_boxes!(ax_left, lit, left)
    _draw_free_pool_points!(ax_right, lit, right)
    _draw_model_marks!(ax_left, lift(p -> p[1:nl], metab_obs))
    _draw_model_marks!(ax_right, lift(p -> p[nl+1:end], metab_obs))
    _metab_legend!(gl[1, 3])
    return (; ax_left, ax_right)
end
```

If Task 2 Step 5 found live cells below 1e-4 µM, change the `_METAB_YLIMS` docstring to state that.

- [ ] **Step 4: Wire it into the explorer**

In `src/explorer.jl`, replace lines 74-83 (from the `# Footer row 2, ...` comment to the end of the `deregister_interaction!` loop) with:

```julia
    # Footer row 2, ATP slider row 3, full-width metabolite panel row 4: two linked axes
    # (free ≈ total | free pool) plus their legend, in a nested GridLayout.
    metab_gl = GridLayout(fig[4, 1:3])
    rowsize!(fig.layout, 4, Relative(0.30))
    lit = load_metab_lit(metab_lit_path())
    metab_obs = Observable(fill(NaN, length(_PANEL_METAB_ORDER)))
    metab_axes = draw_metab_split!(metab_gl, lit, metab_obs)
    for ax in metab_axes, interaction in (:rectanglezoom, :dragpan, :scrollzoom)
        deregister_interaction!(ax, interaction)
    end
```

- [ ] **Step 5: Run the full suite and look at the figure**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=. scripts/render_readme_screenshot.jl
```

Expected: PASS. Then Read `docs/src/assets/explorer.png` and check:
- both titles fit, and the legend overlaps no axis;
- the "AMP / no data" label is readable;
- the model diamonds sit within the free-pool data for ADP, NAD and NADPH;
- the NADP diamond is near its one point.

If a title overflows, shorten it per spec §4 ("check the screenshot for overflow") and re-render. Do **not** commit the PNG yet; Task 7 owns it.

- [ ] **Step 6: Commit**

```bash
git branch --show-current
git add src/metab_panel.jl src/explorer.jl test/test_metab_panel.jl test/test_explorer.jl
git commit -F - <<'EOF'
feat(metab-panel): split the explorer's metabolite panel into free ~ total and free-pool axes

Plain English: the explorer's metabolite plot is now two plots side by side, as in the
PPP_Experiments basic panel. The left compares ordinary metabolites with metabolomics; the right
compares the protein-bound cofactors with free-pool measurements only, marked by how each value
was obtained. Model values too small or too large to show appear as labelled arrows at the edge
instead of disappearing.

load_metab_lit returns _LitPoint NamedTuples and validates method/side; draw_metab_split!
replaces draw_metab_panel! (pinned x-limits, linked log y over _METAB_YLIMS, legend outside the
data, fixed plot count across reselection).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 7: Docs, screenshot, version, final verification

**Repo:** PhaseMap, branch `feat/free-pool-panel-split`.

**Files:**
- Modify: `docs/src/tutorial.md:38-43, 53-65`, `docs/src/science.md:35`, `docs/src/recompute.md:36-45, 53-66`
- Replace: `docs/src/assets/explorer.png`
- Modify: `Project.toml:3`
- Modify (local, never committed): `CLAUDE.md`

- [ ] **Step 1: Tutorial**

In `docs/src/tutorial.md`, replace the sentence ending at line 43 ("...PentosePhosphatePathway.jl's default NADPH and R5P sinks.") with:

```markdown
cell nearest **physiological demand**: PentosePhosphatePathway.jl's default NADPH and R5P sinks,
with the glycolytic activator F26BP at 5 µM (PentosePhosphatePathway.jl's own default is 0).
```

Replace the whole `### The metabolite panel` section body (lines 55-65) with:

```markdown
Below the network sits a full-width **metabolite panel** of 20 metabolites, split in two. The
model tracks each metabolite as its **free** (unbound) concentration, while most measurements
extract the cell and count the **total** pool.

- **Left, free ≈ total** (`Glucose, G6P, F6P, F16BP, DHAP, GAP, ATP, PGLn, PGA, Ru5P, R5P, X5P,
  S7P, E4P`): metabolites whose free and total pools are about the same, so each one is
  compared with a static **box-and-whisker of extraction metabolomics** from cultured mammalian
  cells. The boxes include drug-treated and glucose-pulse measurements as well as basal ones.
  PGLn has no measurements, so it has no box.
- **Right, free pool** (`ADP, AMP, NAD, NADH, NADP, NADPH`): cofactors that are mostly
  protein-bound, so they are compared with **free-pool measurements only**, one point per
  value. Filled points are measured; open points are inferred, by us or by the original
  authors. Circles are sensor-based; diamonds are equilibrium-based. Free NADP⁺ has no direct
  measurement: its one point is HeLa free NADPH divided by the HeLa free NADPH:NADP⁺ ratio. The
  other measured ratios come from cell lines with no free-NADPH measurement, so they cannot be
  paired. AMP has no free-pool measurement, so its column reads "no data".

For the selected cell, each metabolite's steady-state pool (µM of cytosolic water, log y-axis)
is a **purple diamond**. A pool outside the axis range is a purple triangle at the bottom or top
edge, labelled with its value. In ATP-collapsed cells many pools do this at once. The model
marks redraw on every click and slider move; the literature marks are static. GAP and DHAP are
on the free aldehyde/ketone basis on both sides.
```

- [ ] **Step 2: Science page**

In `docs/src/science.md`, after "ATP demand anchors on the glycolytic ceiling `2 · HK1_cap`." (line 35), add:

```markdown
All grids are solved with fructose-2,6-bisphosphate (F26BP), the PFK activator, at 5 µM: the
physiological operating point used across the project. Without it, upper glycolysis runs near
its limit at physiological ATP demand, and ADP builds up about 20-fold.
```

- [ ] **Step 3: Recompute page**

In `docs/src/recompute.md`:

- After the `regenerate_grid()` code block, add the paragraph: "`regenerate_grid(; f26bp_uM = 5.0)` sets the F26BP value (µM of cytosolic water); the default is the phase map's physiological 5 µM, not PentosePhosphatePathway.jl's default of 0."
- In the `[grids]` paragraph, after "the grid axes," insert "the F26BP value,".
- In "Refreshing the literature table", replace "free-pool NADP(H)" with "free-pool cofactors (ADP, AMP, NAD⁺, NADH, NADP⁺, NADPH: free-pool measurements and same-cell-line matched pairs only)".
- Replace the final sentence ("This rewrites ... derived free-NADP⁺ values).") with: "This rewrites `data/literature_pools.csv` and the `[literature]` table (PPP_Experiments commit and source files). Each row's `method` column (`metabolomics`, `sensor` or `equilibrium`) sets its marker in the panel. The export refuses to run if this package's panel grouping (`_FREE_EQ_TOTAL_PANEL`, `_FREE_POOL_PANEL`) differs from PPP_Experiments' `metab_data.jl`."

- [ ] **Step 4: Version, screenshot and local CLAUDE.md**

- `Project.toml` line 3: `version = "0.3.0"`.
- Run `julia --project=. scripts/render_readme_screenshot.jl`, then Read `docs/src/assets/explorer.png` to confirm the split panel.
- Local `CLAUDE.md` (gitignored; do not stage it):
  - In the `grid_io.jl` bullets, change the `row_to_pools` line to "ordered by `_PANEL_METAB_ORDER` = `_FREE_EQ_TOTAL_PANEL` then `_FREE_POOL_PANEL` (grouping guarded by PPP_Experiments' `check_panel_grouping`)". Add "`load_grid` checks renderer columns; `build_explorer` checks pool columns and `atpase_frac`".
  - Replace the `metab_panel.jl` bullet with "`metab_panel.jl` — `draw_metab_split!`: free ≈ total axis (metabolomics boxes) + free-pool axis (marker-keyed points), purple model diamonds with labelled off-scale triangles. The literature table arrives already on the model basis (no corrections here)."
  - Add to the extension notes: "`f26bp_uM` keyword (default `_PHYSIOLOGICAL_DEMAND.f26bp_uM` = 5 µM, recorded in `[grids]`)".

- [ ] **Step 5: Final verification**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=/home/james/projects/PPP_Experiments test/manual/test_extension.jl
julia --project=docs docs/make.jl
git status --short
```

Expected: all three PASS (the docs build has no errors). `git status` shows only the four docs files, `Project.toml` and `explorer.png` (not `CLAUDE.md`).

- [ ] **Step 6: Commit**

```bash
git branch --show-current
git add docs/src/tutorial.md docs/src/science.md docs/src/recompute.md docs/src/assets/explorer.png Project.toml
git commit -F - <<'EOF'
docs: describe the split metabolite panel and the 5 uM F26BP operating point; bump to 0.3.0

Plain English: the tutorial, science and recompute pages now explain the two-part metabolite
panel, the marker key, why NADP has one point and AMP none, and that the phase map runs with
F26BP at 5 uM. The README screenshot shows the new panel. Version 0.3.0 because the grid and
literature file formats changed and old grids no longer open in the explorer.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

- [ ] **Step 7: Finish the branch**

Use superpowers:finishing-a-development-branch. Opening the PhaseMap PR is outward-facing, so ask the user first. The PR description carries the Task 2 F26BP summary table and links the PPP_Experiments PR.
