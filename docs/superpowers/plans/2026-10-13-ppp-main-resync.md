# Resync to PentosePhosphatePathway.jl `main` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the recompute extension's initial state, regenerate both shipped grids from PPP.jl `main` (`fcccae8`), put the metabolite panel's literature boxes on the model's pool basis via a table exported from PPP_Experiments, and record provenance for all shipped data.

**Architecture:** Two repos. PhaseMap (`/home/james/projects/GlycolysisPentosePathwayPhaseMap.jl`, branch `feat/ppp-main-resync`) gains a tiny core provenance module (TOML read-modify-write + git state), a fixed extension, a long-format literature loader, and regenerated data. PPP_Experiments (`/home/james/projects/PPP_Experiments`, branch `feat/phasemap-literature-export`) gains a pure table builder plus a driver that reuses `import_data.jl`/`metab_data.jl` (the only place basis conversion lives) and writes `literature_pools.csv` + the `[literature]` provenance table into PhaseMap's `data/`. All model runs use the PPP_Experiments environment, whose `[sources]` already `develop` both `../GlycolysisPentosePathwayPhaseMap.jl` and `../PentosePhosphatePathway.jl` by path (this is the spec §1 "scratch environment that develops both").

**Tech Stack:** Julia 1.12 (package compat 1.10), Makie/CairoMakie/GLMakie, CSV.jl, DataFrames.jl, TOML stdlib, Distributions (PPP_Experiments only), OrdinaryDiffEq `Rodas5P`.

**Spec:** `docs/superpowers/specs/2026-10-13-ppp-main-resync-design.md` (PhaseMap repo).

## Global Constraints

- PPP.jl checkout `/home/james/projects/PentosePhosphatePathway.jl` must be on `main` at `fcccae8a196ce9cfba23fae91ca77c4be4540807`, clean, before any model run; it stays on `main` afterwards.
- Solver settings unchanged: `Rodas5P()`, `TerminateSteadyState(1e-13, 1e-8)`, `abstol = 1e-14`, `reltol = 1e-9`.
- Grid axes unchanged: `n_nadph = 12`, `n_r5p = 12`, `atp_fracs = 10 .^ range(log10(0.01), log10(0.20), length = 11)`, `atpase_frac = 0.10`, `phi_lo = 1e-3`, `phi_hi = 1.0`.
- `literature_pools.csv` columns exactly `metabolite, value_uM, reference, cell_type, derived`; ASCII only; µM on the model basis (cytosolic water; free NADP(H)).
- All literature rows are kept, including drug-treated and glucose-pulse conditions.
- NADP rows: `nadp_draws = 1000`, `nadp_seed = 20261013`, reference `"derived: free NADPH / R (docs/free_vs_total_measurements.md sec. 4)"`, cell type `"generic mammalian (derived)"`, `derived = true`.
- `data/provenance.toml`: `[grids]` written only by `regenerate_grid`, `[literature]` only by the PPP_Experiments export; each writer replaces only its own table.
- PhaseMap version 0.1.0 → 0.2.0. New dep `TOML` with compat `"1"`.
- 6-phosphogluconate is `PGA` everywhere; never `6PG`.
- Section-header banners: plain ASCII `#`, 90 characters wide, title centred (match `ext/PentosePhosphatePathwayExt.jl`).
- PPP_Experiments: stage named files only (never `git add -A`); confirm `git branch --show-current` before every commit. Commit subjects `type(scope): summary`; bodies open with a plain-English summary; end with the attribution lines below.
- Every commit message ends with:
  ```
  Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01KBx8vahm6QPQeQjsQuYbtw
  ```
- GLMakie needs an X display even to load. If a command that loads PhaseMap fails with a GLFW/display error, prefix it with `xvfb-run -a`.
- Local `CLAUDE.md` in PhaseMap is gitignored: edit it, never commit it.

## Review Focus

1. **Grids silently built from the wrong PPP.jl commit** (someone switches the shared checkout to a test branch). Expected: provenance records commit, branch and dirty flag, and the full run is gated on `fcccae8`/clean. Pinned by Task 2's manual extension test (provenance commit == `git rev-parse HEAD` of `pkgdir(PPP)`) and Task 5 Step 1's pre-run gate + Step 4's post-run assert.
2. **A second provenance writer clobbering the other's table** (regenerating grids after the literature export, or vice versa). Expected: the other table survives byte-for-byte in content. Pinned by Task 1's `update_provenance!` preservation test.
3. **Non-ASCII or invalid bytes in the exported CSV** (Latin-1 source CSV, `§`/`µ` in reference strings). Expected: ASCII-only file. Pinned by Task 3's `isascii` test on a written table and Task 4's shipped-file `isascii` test.
4. **A literature value that is zero, negative or non-finite after conversion** (would break the log-y panel or silently vanish). Expected: the export refuses to write. Pinned by Task 3's `ArgumentError` test on a `0.0` value and Task 4's shipped-table positivity test.
5. **Free-pool data gaining a direct NADP⁺ row, or losing NADPH** (would make the derived-only NADP rule wrong or empty). Expected: the export errors with a message naming the rule. Pinned by Task 3's two guard tests.

---

## File map

PhaseMap (`/home/james/projects/GlycolysisPentosePathwayPhaseMap.jl`):
- Create `src/provenance.jl` — `provenance_path`, `read_provenance`, `update_provenance!`, `git_state`.
- Modify `src/GlycolysisPentosePathwayPhaseMap.jl` — `import TOML`, include `provenance.jl`, `metab_lit_path()` → `literature_pools.csv`.
- Modify `Project.toml` — `TOML` dep + compat, version 0.2.0.
- Modify `ext/PentosePhosphatePathwayExt.jl` — init fix, `[grids]` provenance.
- Modify `src/metab_panel.jl` — long-format loader, drop GAP/DHAP constants, `_METAB_YLABEL`.
- Create `test/test_provenance.jl`; modify `test/runtests.jl`, `test/test_metab_panel.jl`, `test/test_explorer.jl`, `test/test_grid_io.jl`.
- Create `test/manual/test_extension.jl` (run by hand from the PPP_Experiments env; not in `runtests.jl`).
- Data: create `data/literature_pools.csv`, `data/provenance.toml`; delete `data/Metabolite_Concentrations.csv`; regenerate `data/atp_demand_grid.csv`, `data/demand_grid.csv`.
- Docs: `docs/src/recompute.md`, `docs/src/tutorial.md`, `docs/src/index.md`, `docs/src/assets/explorer.png`; local `CLAUDE.md`.

PPP_Experiments (`/home/james/projects/PPP_Experiments`):
- Create `scripts/phasemap_export/literature_export.jl` — pure `literature_pools_table`.
- Create `scripts/phasemap_export/test_literature_export.jl`.
- Create `scripts/phasemap_export/export_literature_pools.jl` — driver.
- Modify `scripts/README.md` — one entry for `phasemap_export/`.

---

### Task 1: Core provenance helpers (PhaseMap)

**Files:**
- Create: `src/provenance.jl`
- Modify: `src/GlycolysisPentosePathwayPhaseMap.jl`, `Project.toml`
- Test: `test/test_provenance.jl`, `test/runtests.jl`

**Interfaces:**
- Consumes: `DATA_DIR` (already in the main module).
- Produces (all unexported, used as `GlycolysisPentosePathwayPhaseMap.<name>`):
  - `provenance_path() -> String` = `joinpath(DATA_DIR, "provenance.toml")`
  - `read_provenance(path = provenance_path()) -> Dict{String,Any}`
  - `update_provenance!(path::AbstractString, section::AbstractString, table::AbstractDict) -> Dict{String,Any}` — replaces only `section`, creates the file if absent, writes sorted keys.
  - `git_state(dir::AbstractString; paths = String[]) -> Union{Nothing, NamedTuple{(:commit, :branch, :dirty), Tuple{String,String,Bool}}}` — `nothing` unless `dir` is the top level of a git work tree; `dirty` restricted to `paths` when non-empty.

- [x] **Step 1: Add the TOML dependency**

In `Project.toml`, add under `[deps]` (keep alphabetical order):
```toml
TOML = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
```
and under `[compat]`:
```toml
TOML = "1"
```
Then resolve both environments that load this package:
```bash
cd /home/james/projects/GlycolysisPentosePathwayPhaseMap.jl && julia --project=. -e 'using Pkg; Pkg.resolve()'
julia --project=/home/james/projects/PPP_Experiments -e 'using Pkg; Pkg.resolve()'
```
Expected: both finish without error (Manifests are gitignored in both repos).

- [x] **Step 2: Write the failing test**

Create `test/test_provenance.jl`:
```julia
using Test, GlycolysisPentosePathwayPhaseMap
const G = GlycolysisPentosePathwayPhaseMap

@testset "update_provenance! creates, replaces one table, preserves the other" begin
    path = joinpath(mktempdir(), "provenance.toml")
    G.update_provenance!(path, "grids", Dict("ppp_commit" => "abc", "atp_grid_rows" => 3))
    @test G.read_provenance(path)["grids"]["ppp_commit"] == "abc"

    G.update_provenance!(path, "literature", Dict("ppp_experiments_commit" => "def"))
    p = G.read_provenance(path)
    @test p["grids"] == Dict("ppp_commit" => "abc", "atp_grid_rows" => 3)
    @test p["literature"]["ppp_experiments_commit"] == "def"

    G.update_provenance!(path, "grids", Dict("ppp_commit" => "xyz"))
    p = G.read_provenance(path)
    @test p["grids"] == Dict("ppp_commit" => "xyz")          # old keys dropped
    @test p["literature"] == Dict("ppp_experiments_commit" => "def")   # untouched
end

@testset "git_state" begin
    @test G.git_state(mktempdir()) === nothing               # not a git work tree
    repo = mktempdir()
    run(pipeline(`git -C $repo init -q`; stderr = devnull))
    write(joinpath(repo, "a.txt"), "1")
    run(`git -C $repo -c user.name=t -c user.email=t@t add a.txt`)
    run(`git -C $repo -c user.name=t -c user.email=t@t commit -q -m init`)
    s = G.git_state(repo)
    @test length(s.commit) == 40
    @test !s.dirty
    write(joinpath(repo, "b.txt"), "2")
    @test G.git_state(repo).dirty                            # untracked file -> dirty
    @test !G.git_state(repo; paths = ["a.txt"]).dirty        # scoped to clean path
    @test G.git_state(joinpath(repo, "sub_does_not_exist")) === nothing
end
```
Add `include("test_provenance.jl")` to `test/runtests.jl` as the first include inside the testset.

- [x] **Step 3: Run test to verify it fails**

Run: `cd /home/james/projects/GlycolysisPentosePathwayPhaseMap.jl && julia --project=. -e 'using GlycolysisPentosePathwayPhaseMap, Test; include("test/test_provenance.jl")'`
Expected: FAIL/ERROR with `UndefVarError: update_provenance! not defined`.

- [x] **Step 4: Implement**

Create `src/provenance.jl`:
```julia
# provenance.jl — where the shipped data came from (data/provenance.toml)
#
# Two writers share one file, each owning one table: `regenerate_grid` (the recompute
# extension) writes [grids]; the PPP_Experiments literature export writes [literature].
# update_provenance! read-modify-writes only the named table so neither clobbers the other.

"Path of the shipped provenance file."
provenance_path() = joinpath(DATA_DIR, "provenance.toml")

"Parse a provenance TOML file into a `Dict{String,Any}`."
read_provenance(path = provenance_path()) = TOML.parsefile(path)

"""
    update_provenance!(path, section, table) -> Dict{String,Any}

Replace table `section` of the provenance TOML at `path` with `table` (creating the file if
absent), leaving every other table untouched. Keys are written sorted. Returns the full
parsed content as written.
"""
function update_provenance!(path::AbstractString, section::AbstractString, table::AbstractDict)
    prov = isfile(path) ? TOML.parsefile(path) : Dict{String,Any}()
    prov[String(section)] = Dict{String,Any}(String(k) => v for (k, v) in table)
    mkpath(dirname(path))
    open(io -> TOML.print(io, prov; sorted = true), path, "w")
    return prov
end

"""
    git_state(dir; paths = String[]) -> Union{Nothing, NamedTuple}

`(; commit, branch, dirty)` of the git work tree whose top level is `dir`, or `nothing` when
`dir` is not a work-tree top level (e.g. a registry-installed package). `dirty` counts
untracked files too and is restricted to `paths` (relative to `dir`) when given.
"""
function git_state(dir::AbstractString; paths = String[])
    isdir(dir) || return nothing
    git(args...) = readchomp(pipeline(`git -C $dir $args`; stderr = devnull))
    top = try git("rev-parse", "--show-toplevel") catch; return nothing end
    realpath(top) == realpath(dir) || return nothing
    commit = git("rev-parse", "HEAD")
    branch = git("rev-parse", "--abbrev-ref", "HEAD")
    dirty  = !isempty(git("status", "--porcelain", "--", paths...))
    return (; commit, branch, dirty)
end
```
In `src/GlycolysisPentosePathwayPhaseMap.jl`, change the imports and includes to:
```julia
using Makie, CSV, DataFrames
import GLMakie, CairoMakie   # backends activated per-entry; NOT `using` (avoids dual re-export clash)
import TOML
```
and add `include("provenance.jl")` immediately before `include("grid_io.jl")`.

- [x] **Step 5: Run test to verify it passes**

Run: `cd /home/james/projects/GlycolysisPentosePathwayPhaseMap.jl && julia --project=. -e 'using GlycolysisPentosePathwayPhaseMap, Test; include("test/test_provenance.jl")'`
Expected: both testsets PASS.

- [x] **Step 6: Commit**

```bash
cd /home/james/projects/GlycolysisPentosePathwayPhaseMap.jl && git branch --show-current   # must print feat/ppp-main-resync
git add Project.toml src/provenance.jl src/GlycolysisPentosePathwayPhaseMap.jl test/test_provenance.jl test/runtests.jl
git commit -F - <<'EOF'
feat(provenance): add a shared provenance file with per-writer tables

Adds a small helper that records where the package's shipped data came from
(which model commit built the grids, which data commit built the literature
table) in data/provenance.toml, with each producer updating only its own part.

- src/provenance.jl: provenance_path, read_provenance, update_provenance!
  (read-modify-write of one TOML table), git_state (commit/branch/dirty of a
  work-tree top level, optionally path-scoped).
- TOML stdlib dep (compat "1").

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KBx8vahm6QPQeQjsQuYbtw
EOF
```

---

### Task 2: Extension fix + `[grids]` provenance (PhaseMap)

**Files:**
- Modify: `ext/PentosePhosphatePathwayExt.jl` (`_demand_setup` line `init0 = …`; `_regenerate_grid_impl`)
- Test: `test/manual/test_extension.jl` (manual; run from the PPP_Experiments env)

**Interfaces:**
- Consumes: `GPPM.git_state`, `GPPM.update_provenance!`, `GPPM.read_provenance` (Task 1); `PPP.glycolysis_ppp_init_conc`.
- Produces: `regenerate_grid(; atp_out, demand_out, provenance_out = joinpath(dirname(atp_out), "provenance.toml"), n_nadph, n_r5p, atp_fracs, atpase_frac, phi_lo, phi_hi) -> (; atp_out, demand_out, provenance_out)`; `[grids]` keys: `ppp_commit`, `ppp_branch`, `ppp_dirty` (git checkouts only), `cellmetabolism_version`, `generated`, `n_nadph`, `n_r5p`, `atp_fracs`, `atpase_frac`, `phi_lo`, `phi_hi`, `atp_grid_rows`, `demand_grid_rows`, `atp_grid_nonconverged`, `demand_grid_nonconverged`.

- [x] **Step 1: Gate on the PPP.jl checkout**

Run:
```bash
git -C /home/james/projects/PentosePhosphatePathway.jl rev-parse HEAD
git -C /home/james/projects/PentosePhosphatePathway.jl status --porcelain
git -C /home/james/projects/PentosePhosphatePathway.jl branch --show-current
```
Expected: `fcccae8a196ce9cfba23fae91ca77c4be4540807`, empty output, `main`. If not, STOP and report.

- [x] **Step 2: Write the failing manual test**

Create `test/manual/test_extension.jl`:
```julia
# Manual extension check. PPP.jl is unregistered, so CI cannot load the extension; run this by
# hand from the PPP_Experiments environment, which `develop`s this package and PPP.jl by path:
#   julia --project=/home/james/projects/PPP_Experiments test/manual/test_extension.jl
using Test, CSV, DataFrames
using GlycolysisPentosePathwayPhaseMap, PentosePhosphatePathway
const G   = GlycolysisPentosePathwayPhaseMap
const PPP = PentosePhosphatePathway

ext = Base.get_extension(G, :PentosePhosphatePathwayExt)

@testset "recompute extension" begin
    @test ext !== nothing

    # initial state = PPP.jl's canonical combined init (with PPP_glycolysis_init_overrides)
    s = ext._demand_setup(0.10, 2, 2, 1e-3, 1.0)
    @test keys(s.init0) == keys(PPP.glycolysis_ppp_init_conc)
    @test collect(s.init0) == collect(PPP.glycolysis_ppp_init_conc)
    @test s.init0 !== PPP.glycolysis_ppp_init_conc          # a copy, never the package global

    # 2 x 2 smoke regenerate into scratch paths, with provenance
    dir = mktempdir()
    out = regenerate_grid(; atp_out = joinpath(dir, "atp.csv"),
                            demand_out = joinpath(dir, "demand.csv"),
                            n_nadph = 2, n_r5p = 2, atp_fracs = [0.10])
    @test out.provenance_out == joinpath(dir, "provenance.toml")
    atp = G.load_grid(out.atp_out)
    @test nrow(atp) == 4
    @test count(==(:Terminated), atp.retcode) >= 1
    prov = G.read_provenance(out.provenance_out)["grids"]
    gs = G.git_state(pkgdir(PPP))
    @test prov["ppp_commit"] == gs.commit
    @test prov["ppp_branch"] == gs.branch
    @test prov["ppp_dirty"] == gs.dirty
    @test prov["atp_grid_rows"] == 4
    @test prov["demand_grid_rows"] == 4
    @test prov["atp_grid_nonconverged"] == count(!=(:Terminated), atp.retcode)
    @test prov["atp_fracs"] == [0.10]
    @test !haskey(G.read_provenance(out.provenance_out), "literature")
end
```

- [x] **Step 3: Run it to verify it fails**

Run: `cd /home/james/projects/GlycolysisPentosePathwayPhaseMap.jl && julia --project=/home/james/projects/PPP_Experiments test/manual/test_extension.jl`
Expected: FAIL on `keys(s.init0) == keys(PPP.glycolysis_ppp_init_conc)` or `collect(...) ==` (hand-merge skips the overrides), and an error on `out.provenance_out` (field missing). (First run precompiles for several minutes; that is not a hang.)

- [x] **Step 4: Implement**

In `ext/PentosePhosphatePathwayExt.jl`, `_demand_setup`, replace
```julia
    init0   = PPP.merge_LArrays(PPP.PPP_init_conc, CellMetabolism.glycolysis_init_conc)
```
with
```julia
    # PPP.jl's canonical combined init: PPP_init_conc, then PPP_glycolysis_init_overrides
    # (FIX-15/16/22), then CellMetabolism's. A hand-merge skips the overrides. Copied: the
    # export is a package-level array that the warm-start sweep must never write into.
    init0   = copy(PPP.glycolysis_ppp_init_conc)
```
Replace the whole `_regenerate_grid_impl` function with:
```julia
# Real implementation behind GPPM.regenerate_grid. Registered into the core module's hook Ref by
# __init__ (below) at load time — never defined as a method on GPPM.regenerate_grid, which would
# be a precompile-breaking method overwrite. Docstring lives on the exported parent function.
function _regenerate_grid_impl(; atp_out = GPPM.default_atp_grid(),
                                 demand_out = GPPM.default_demand_grid(),
                                 provenance_out = joinpath(dirname(atp_out), "provenance.toml"),
                                 n_nadph = 12, n_r5p = 12,
                                 atp_fracs = 10 .^ range(log10(0.01), log10(0.20), length = 11),
                                 atpase_frac = 0.10, phi_lo = 1e-3, phi_hi = 1.0)
    mkpath(dirname(atp_out)); mkpath(dirname(demand_out))
    atp = run_atp_grid(; n_nadph, n_r5p, atp_fracs, phi_lo, phi_hi)
    CSV.write(atp_out, atp)
    demand = run_demand_grid(; n_nadph, n_r5p, atpase_frac, phi_lo, phi_hi)
    CSV.write(demand_out, demand)
    GPPM.update_provenance!(provenance_out, "grids",
        _grid_provenance(atp, demand; n_nadph, n_r5p, atp_fracs, atpase_frac, phi_lo, phi_hi))
    @info "Regenerated grids" atp_out demand_out provenance_out atp_rows=nrow(atp) demand_rows=nrow(demand)
    return (; atp_out, demand_out, provenance_out)
end

# The [grids] provenance table: which model build produced the grids, with what axes, and how
# many cells failed to reach steady state. A git checkout of PPP.jl records commit/branch/dirty;
# a registry install records its version instead.
function _grid_provenance(atp, demand; n_nadph, n_r5p, atp_fracs, atpase_frac, phi_lo, phi_hi)
    gs = GPPM.git_state(pkgdir(PPP))
    prov = Dict{String,Any}(
        "ppp_commit"               => gs === nothing ? "v" * string(pkgversion(PPP)) : gs.commit,
        "cellmetabolism_version"   => string(pkgversion(CellMetabolism)),
        "generated"                => Libc.strftime("%Y-%m-%d", time()),
        "n_nadph" => n_nadph, "n_r5p" => n_r5p,
        "atp_fracs" => collect(Float64, atp_fracs), "atpase_frac" => Float64(atpase_frac),
        "phi_lo" => Float64(phi_lo), "phi_hi" => Float64(phi_hi),
        "atp_grid_rows"            => nrow(atp),
        "demand_grid_rows"         => nrow(demand),
        "atp_grid_nonconverged"    => count(!=(:Terminated), atp.retcode),
        "demand_grid_nonconverged" => count(!=(:Terminated), demand.retcode),
    )
    if gs !== nothing
        prov["ppp_branch"] = gs.branch
        prov["ppp_dirty"]  = gs.dirty
    end
    return prov
end
```
Also update the `regenerate_grid` docstring in `src/GlycolysisPentosePathwayPhaseMap.jl` to:
```julia
"""
    regenerate_grid(; kwargs...)

Recompute the shipped grid CSVs by solving the kinetic ODE model, and record the model build
that produced them in the `[grids]` table of `data/provenance.toml` (next to `atp_out` unless
`provenance_out` is given). Only available when the recompute extension is active — load it
first with `using PentosePhosphatePathway`.
"""
```

- [x] **Step 5: Run it to verify it passes**

Run: `cd /home/james/projects/GlycolysisPentosePathwayPhaseMap.jl && julia --project=/home/james/projects/PPP_Experiments test/manual/test_extension.jl`
Expected: all tests PASS. Also run the core suite: `julia --project=. -e 'using Pkg; Pkg.test()'` → PASS (the core never loads the extension).

- [x] **Step 6: Commit**

```bash
cd /home/james/projects/GlycolysisPentosePathwayPhaseMap.jl && git branch --show-current   # feat/ppp-main-resync
git add ext/PentosePhosphatePathwayExt.jl src/GlycolysisPentosePathwayPhaseMap.jl test/manual/test_extension.jl
git commit -F - <<'EOF'
fix(ext): start sweeps from PPP.jl's canonical initial state and record grid provenance

The grid recompute was building the model's starting concentrations by hand,
which skipped PPP.jl's corrected values for adenylates, phosphate, pyruvate,
DHAP, F16BP and PEP. It now uses the starting state PPP.jl itself publishes,
and every regeneration records which model commit produced the grids.

- _demand_setup: init0 = copy(PPP.glycolysis_ppp_init_conc) (was a
  merge_LArrays hand-merge that skipped PPP_glycolysis_init_overrides).
- _regenerate_grid_impl: new provenance_out kwarg; writes the [grids] table
  (PPP.jl commit/branch/dirty, CellMetabolism version, axes, row and
  non-converged counts).
- test/manual/test_extension.jl: manual check run from the PPP_Experiments
  env (PPP.jl is unregistered, so CI cannot load the extension).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KBx8vahm6QPQeQjsQuYbtw
EOF
```

---

### Task 3: Literature table builder (PPP_Experiments)

**Files:**
- Create: `scripts/phasemap_export/literature_export.jl`
- Test: `scripts/phasemap_export/test_literature_export.jl`

**Interfaces:**
- Consumes (from `scripts/helper_files/metab_data.jl`, included by the caller): `metab_columns(df)`, `FREE_POOL_METABS`, `nadph_nadp_ratio_inputs(free_df)`, `draw_free_nadph_nadp_ratio(rng, inputs)`. Inputs are in **M**: `metab_df` after `import_data.jl` (÷1e6, basis-converted, GAP/DHAP-corrected); `free_df` after `prepare_free_pool` (concentration rows have `Unit == "M"`).
- Produces: `literature_pools_table(metab_df, free_df, panel; nadp_draws = 1000, nadp_seed = 20261013) -> DataFrame` with columns `metabolite::String, value_uM::Float64, reference::String, cell_type::String, derived::Bool`; constants `NADP_DERIVED_REFERENCE`, `NADP_DERIVED_CELL_TYPE`.

- [ ] **Step 1: Write the failing test**

Create `scripts/phasemap_export/test_literature_export.jl`:
```julia
# Run from the repo root: julia --project scripts/phasemap_export/test_literature_export.jl
using Test, CSV, DataFrames, Random, Statistics
include(joinpath(@__DIR__, "..", "helper_files", "metab_data.jl"))
include(joinpath(@__DIR__, "literature_export.jl"))

# metab_df as import_data.jl leaves it: values in M on the model basis.
metab_df() = DataFrame("Reference" => ["Ref A", "Ref B", "Ref C"],
                       "Cell Type" => ["HeLa", "HeLa_Oligomycin", "iBMK"],
                       "Basis" => ["whole_cell", "whole_cell", "cell_water"],
                       "Notes" => ["", "", ""],
                       "G6P" => Union{Missing,Float64}[100e-6, 50e-6, missing],
                       "PGLn" => Union{Missing,Float64}[missing, missing, missing])

# free_df as prepare_free_pool leaves it: concentrations in M; the 5 sec. 4 ratio sources.
function free_df(; extra = DataFrame())
    df = DataFrame("Reference" => ["Tao", "Moon", "Mito", "Sallin 1", "Sallin 2", "Sallin 3",
                                   "Sallin 4", "Veech"],
                   "Cell Type" => ["HeLa", "HeLa", "HeLa", "HEK293T", "NIH/3T3", "U2OS",
                                   "HeLa", "rat liver"],
                   "Compartment" => ["cytosol", "cytosol", "mitochondria", "cytosol",
                                     "cytosol", "cytosol", "cytosol", "cytosol"],
                   "Quantity" => ["NADPH", "NADPH", "NADPH", "NADPH/NADP", "NADPH/NADP",
                                  "NADPH/NADP", "NADPH/NADP", "NADPH/NADP"],
                   "Value" => [3.1e-6, 2.6e-6, 9.9e-6, 21.6, 39.5, 55.8, 75.0, 91.6],
                   "Unit" => ["M", "M", "M", "ratio", "ratio", "ratio", "ratio", "ratio"],
                   "Method" => ["sensor", "sensor", "sensor", "sensor", "sensor", "sensor",
                                "sensor", "equilibrium"])
    return vcat(df, extra; cols = :union)
end

const PANEL = (:G6P, :PGLn, :NADPH, :NADP)

@testset "literature_pools_table" begin
    t = literature_pools_table(metab_df(), free_df(), PANEL; nadp_draws = 200)
    @test names(t) == ["metabolite", "value_uM", "reference", "cell_type", "derived"]

    g = t[t.metabolite .== "G6P", :]
    @test g.value_uM ≈ [100.0, 50.0]                          # M -> µM; missing dropped
    @test g.cell_type == ["HeLa", "HeLa_Oligomycin"]          # treated rows kept
    @test !any(g.derived)
    @test !any(t.metabolite .== "PGLn")                       # no data -> no rows

    n = t[t.metabolite .== "NADPH", :]
    @test n.value_uM ≈ [3.1, 2.6]                             # cytosol only (mito dropped)
    @test n.reference == ["Tao", "Moon"]

    d = t[t.metabolite .== "NADP", :]
    @test nrow(d) == 200 && all(d.derived)
    @test all(==(NADP_DERIVED_REFERENCE), d.reference)
    @test all(x -> isfinite(x) && x > 0, d.value_uM)
    @test 0.005 < median(d.value_uM) < 0.5                    # ~3 µM / R(~50) ~ 0.06 µM

    t2 = literature_pools_table(metab_df(), free_df(), PANEL; nadp_draws = 200)
    @test t2.value_uM == t.value_uM                           # seeded -> reproducible
    t3 = literature_pools_table(metab_df(), free_df(), PANEL; nadp_draws = 200, nadp_seed = 1)
    @test t3.value_uM[t3.metabolite .== "NADP"] != d.value_uM

    path = tempname() * ".csv"
    CSV.write(path, t)
    @test isascii(read(path, String))
end

@testset "literature_pools_table guards" begin
    @test_throws ArgumentError literature_pools_table(metab_df(), free_df(), (:Foo,))
    bad = metab_df(); bad.G6P[1] = 0.0
    @test_throws ArgumentError literature_pools_table(bad, free_df(), (:G6P,))
    direct = DataFrame("Reference" => ["X"], "Cell Type" => ["HeLa"], "Compartment" => ["cytosol"],
                       "Quantity" => ["NADP"], "Value" => [1e-8], "Unit" => ["M"],
                       "Method" => ["sensor"])
    @test_throws ArgumentError literature_pools_table(metab_df(), free_df(; extra = direct), (:NADP,))
    no_nadph = filter(r -> r.Quantity != "NADPH", free_df())
    @test_throws ArgumentError literature_pools_table(metab_df(), no_nadph, (:NADP,))
    nonascii = metab_df(); nonascii.Reference[1] = "Müller"
    @test_throws ArgumentError literature_pools_table(nonascii, free_df(), (:G6P,))
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /home/james/projects/PPP_Experiments && julia --project scripts/phasemap_export/test_literature_export.jl`
Expected: ERROR `could not open file .../literature_export.jl`.

- [ ] **Step 3: Implement**

Create `scripts/phasemap_export/literature_export.jl`:
```julia
##########################################################################################
#           literature_export.jl -- PhaseMap metabolite-panel literature table           #
##########################################################################################
# Builds the long-format table GlycolysisPentosePathwayPhaseMap.jl ships as
# data/literature_pools.csv: one row per literature value, in µM on the model's pool basis.
# Inputs are the frames import_data.jl leaves behind (values in M, already converted by row
# Basis and species class, GAP/DHAP free-aldehyde corrected), so the basis conversion lives
# only in metab_data.jl. Caller must have included metab_data.jl.
#
# - FREE_POOL_METABS on the panel (NADPH, NADP) come from the free-pool table, never from the
#   total-pool metabolomics columns.
# - NADP has no direct free measurement: it is free NADPH / R, R drawn from the
#   docs/free_vs_total_measurements.md sec. 4 sampler, seeded so the file is reproducible.
# - Every metabolomics row is kept, treated conditions (<line>_<treatment>) included.

const NADP_DERIVED_REFERENCE = "derived: free NADPH / R (docs/free_vs_total_measurements.md sec. 4)"
const NADP_DERIVED_CELL_TYPE = "generic mammalian (derived)"

"Cytosolic free-pool concentration rows (M) of metabolite `m`."
_free_rows(free_df, m::Symbol) = free_df[(free_df.Compartment .== "cytosol") .&
                                         (free_df.Unit .== "M") .&
                                         (free_df.Quantity .== string(m)), :]

"""
    literature_pools_table(metab_df, free_df, panel; nadp_draws = 1000, nadp_seed = 20261013)
        -> DataFrame

Rows `metabolite, value_uM, reference, cell_type, derived` for every metabolite in `panel`.
Throws `ArgumentError` when a panel metabolite has no source (no metabolomics column and not a
free-pool metabolite), when direct free NADP rows exist (the derived-only NADP rule would be
stale), when NADP is requested but there is no free NADPH, when any value is not finite and
> 0, or when any text field is not ASCII.
"""
function literature_pools_table(metab_df, free_df, panel; nadp_draws::Int = 1000,
                                nadp_seed::Integer = 20261013)
    out = DataFrame(metabolite = String[], value_uM = Float64[], reference = String[],
                    cell_type = String[], derived = Bool[])
    mcols = metab_columns(metab_df)
    for m in panel
        if m === :NADP
            isempty(_free_rows(free_df, :NADP)) || throw(ArgumentError(
                "free-pool data now has direct free NADP rows; revisit the derived-only NADP " *
                "rule in literature_export.jl before exporting"))
            nadph = _free_rows(free_df, :NADPH).Value
            isempty(nadph) && throw(ArgumentError(
                "NADP is derived as free NADPH / R, but the free-pool data has no cytosolic NADPH"))
            rng = Random.Xoshiro(nadp_seed)
            inputs = nadph_nadp_ratio_inputs(free_df)
            for _ in 1:nadp_draws
                v = rand(rng, nadph) / draw_free_nadph_nadp_ratio(rng, inputs)
                push!(out, ("NADP", v * 1e6, NADP_DERIVED_REFERENCE, NADP_DERIVED_CELL_TYPE, true))
            end
        elseif m in FREE_POOL_METABS
            for r in eachrow(_free_rows(free_df, m))
                push!(out, (string(m), r.Value * 1e6, String(r.Reference),
                            String(r["Cell Type"]), false))
            end
        elseif string(m) in mcols
            for r in eachrow(metab_df)
                v = r[string(m)]
                ismissing(v) && continue
                push!(out, (string(m), v * 1e6, String(r.Reference), String(r["Cell Type"]), false))
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

- [ ] **Step 4: Run test to verify it passes**

Run: `cd /home/james/projects/PPP_Experiments && julia --project scripts/phasemap_export/test_literature_export.jl`
Expected: both testsets PASS.

- [ ] **Step 5: Commit**

```bash
cd /home/james/projects/PPP_Experiments && git branch --show-current   # must print feat/phasemap-literature-export
git add scripts/phasemap_export/literature_export.jl scripts/phasemap_export/test_literature_export.jl
git commit -F - <<'EOF'
feat(phasemap-export): build the phase-map literature table on the model basis

Adds the function that turns this repo's literature metabolite data into the
table the phase-map package shows behind its model points, already converted
to the model's basis (cytosolic water, free NADP(H)), so the package itself
does no unit or basis conversion.

- literature_pools_table(metab_df, free_df, panel): metabolomics columns for
  free ~ total species, cytosolic free-pool rows for NADPH, and NADP derived
  as free NADPH / R from the free_vs_total sec. 4 sampler (seeded).
- Guards: unknown panel metabolite, direct free NADP rows, missing NADPH,
  non-positive/non-finite values, non-ASCII text.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KBx8vahm6QPQeQjsQuYbtw
EOF
```

---

### Task 4: Export driver + PhaseMap long-format loader and shipped table

**Files:**
- Create (PPP_Experiments): `scripts/phasemap_export/export_literature_pools.jl`; modify `scripts/README.md`
- Modify (PhaseMap): `src/metab_panel.jl`, `src/GlycolysisPentosePathwayPhaseMap.jl` (`metab_lit_path`), `test/test_metab_panel.jl`, `test/test_explorer.jl`
- Data (PhaseMap): create `data/literature_pools.csv`, `data/provenance.toml` (`[literature]`); delete `data/Metabolite_Concentrations.csv`

**Interfaces:**
- Consumes: `literature_pools_table` (Task 3); `GPPM.update_provenance!`, `GPPM.git_state`, `GPPM._PANEL_METAB_ORDER` (PhaseMap core); `import_data.jl` globals `metab_conc_data`, `free_pool_data`.
- Produces: `load_metab_lit(csv_path; order = _PANEL_METAB_ORDER) -> Dict{Symbol, Vector{Float64}}` reading the long format; `_METAB_YLABEL::String`; `[literature]` keys `ppp_experiments_commit`, `ppp_experiments_branch`, `ppp_experiments_dirty`, `source_files`, `nadp_seed`, `nadp_draws`, `generated`, `rows`.

- [ ] **Step 1: Write the failing PhaseMap tests**

Replace the `"load_metab_lit"` testset in `test/test_metab_panel.jl` with (and change the first line to `using Test, CSV, CairoMakie, DataFrames, GlycolysisPentosePathwayPhaseMap`):
```julia
@testset "load_metab_lit (shipped literature_pools.csv)" begin
    path = G.metab_lit_path()
    @test basename(path) == "literature_pools.csv"
    lit = G.load_metab_lit(path)
    @test Set(keys(lit)) == Set(G._PANEL_METAB_ORDER)
    @test all(v -> v isa Vector{Float64}, values(lit))
    @test all(v -> all(x -> isfinite(x) && x > 0, v), values(lit))
    @test all(s -> !isempty(lit[s]), filter(!=(:PGLn), collect(G._PANEL_METAB_ORDER)))
    raw = CSV.read(path, DataFrame)
    @test names(raw) == ["metabolite", "value_uM", "reference", "cell_type", "derived"]
    @test all(raw.derived[raw.metabolite .== "NADP"])
    @test !any(raw.derived[raw.metabolite .!= "NADP"])
    @test isascii(read(path, String))
    @test !isfile(joinpath(G.DATA_DIR, "Metabolite_Concentrations.csv"))
    @test !isdefined(G, :_GAP_CORRECTION) && !isdefined(G, :_DHAP_CORRECTION)
end

@testset "load_metab_lit long format (synthetic)" begin
    path = tempname() * ".csv"
    CSV.write(path, DataFrame(metabolite = ["G6P", "G6P", "Foo"], value_uM = [1.0, 2.0, 3.0],
                              reference = ["a", "b", "c"], cell_type = ["x", "y", "z"],
                              derived = [false, false, false]))
    lit = G.load_metab_lit(path)
    @test lit[:G6P] == [1.0, 2.0]          # no conversion applied by the loader
    @test isempty(lit[:PGLn])
    @test !haskey(lit, :Foo)
end
```
In the `"draw_metab_panel! offscreen"` testset add after `@test ax.yscale[] === log10`:
```julia
    @test ax.ylabel[] == G._METAB_YLABEL
    @test occursin("cytosolic water", G._METAB_YLABEL) && occursin("free", G._METAB_YLABEL)
```
In `test/test_explorer.jl` replace
```julia
    ax_metab = only(filter(x -> x isa Axis && x.ylabel[] == "[metabolite], µM", fig.content))
```
with
```julia
    ax_metab = only(filter(x -> x isa Axis && x.ylabel[] == G._METAB_YLABEL, fig.content))
```

- [ ] **Step 2: Run to verify failure**

Run: `cd /home/james/projects/GlycolysisPentosePathwayPhaseMap.jl && julia --project=. -e 'using GlycolysisPentosePathwayPhaseMap, Test; include("test/test_metab_panel.jl")'`
Expected: FAIL (`basename(path) == "literature_pools.csv"`, `_METAB_YLABEL` undefined).

- [ ] **Step 3: Implement the PhaseMap loader**

In `src/GlycolysisPentosePathwayPhaseMap.jl`:
```julia
metab_lit_path()      = joinpath(DATA_DIR, "literature_pools.csv")
```
In `src/metab_panel.jl`, delete the GAP/DHAP comment block and both constants (lines 3-9), and replace `load_metab_lit` with:
```julia
# The literature table is exported by PPP_Experiments (scripts/phasemap_export/) already on the
# model's pool basis: µM of cytosolic water, free NADP(H) (NADP derived as free NADPH / R),
# GAP/DHAP as the free aldehyde. This package does no basis conversion.

"Metabolite-panel y-axis label; states the basis both the boxes and the diamonds are on."
const _METAB_YLABEL = "[metabolite], µM (cytosolic water; NADP(H) free)"

"""
    load_metab_lit(csv_path; order = _PANEL_METAB_ORDER) -> Dict{Symbol, Vector{Float64}}

Read the long-format literature table (`metabolite, value_uM, reference, cell_type, derived`)
and return, per panel metabolite, its values in µM. Values are used as shipped: the table is
already on the model's pool basis. A panel metabolite with no rows maps to an empty vector;
rows for metabolites outside `order` are ignored.
"""
function load_metab_lit(csv_path; order = _PANEL_METAB_ORDER)
    df = CSV.read(csv_path, DataFrame; types = Dict(:metabolite => String, :value_uM => Float64))
    lit = Dict{Symbol, Vector{Float64}}(sym => Float64[] for sym in order)
    for r in eachrow(df)
        sym = Symbol(r.metabolite)
        haskey(lit, sym) && push!(lit[sym], r.value_uM)
    end
    return lit
end
```
In `draw_metab_panel!`, replace `ax.ylabel = "[metabolite], µM"` with `ax.ylabel = _METAB_YLABEL`.

- [ ] **Step 4: Write the export driver (PPP_Experiments)**

Create `scripts/phasemap_export/export_literature_pools.jl`:
```julia
##########################################################################################
#        export_literature_pools.jl -- refresh PhaseMap's data/literature_pools.csv        #
##########################################################################################
# Run from the PPP_Experiments root (the env `develop`s GlycolysisPentosePathwayPhaseMap.jl):
#   julia --project scripts/phasemap_export/export_literature_pools.jl [<phasemap data dir>]
# Default target: ../GlycolysisPentosePathwayPhaseMap.jl/data. Writes literature_pools.csv and
# the [literature] table of provenance.toml there (the [grids] table is left untouched).

using PentosePhosphatePathway, CSV, DataFrames, Random
using GlycolysisPentosePathwayPhaseMap
const GPPM = GlycolysisPentosePathwayPhaseMap

include(joinpath(@__DIR__, "..", "helper_files", "import_data.jl"))  # metab_conc_data, free_pool_data
include(joinpath(@__DIR__, "literature_export.jl"))

const NADP_DRAWS = 1000
const NADP_SEED  = 20261013
const REPO_ROOT  = normpath(joinpath(@__DIR__, "..", ".."))
const SOURCE_FILES = ["data/Metabolite_Concentrations.csv", "data/Free_Pool_Concentrations.csv",
                      "scripts/helper_files/import_data.jl", "scripts/helper_files/metab_data.jl",
                      "scripts/phasemap_export/literature_export.jl",
                      "scripts/phasemap_export/export_literature_pools.jl"]

out_dir = isempty(ARGS) ? normpath(joinpath(REPO_ROOT, "..", "GlycolysisPentosePathwayPhaseMap.jl", "data")) :
                          ARGS[1]
isdir(out_dir) || error("PhaseMap data dir not found: $out_dir")

tbl = literature_pools_table(metab_conc_data, free_pool_data, GPPM._PANEL_METAB_ORDER;
                             nadp_draws = NADP_DRAWS, nadp_seed = NADP_SEED)
csv_out = joinpath(out_dir, "literature_pools.csv")
CSV.write(csv_out, tbl)

gs = GPPM.git_state(REPO_ROOT; paths = SOURCE_FILES)
gs === nothing && error("PPP_Experiments at $REPO_ROOT is not a git work tree")
GPPM.update_provenance!(joinpath(out_dir, "provenance.toml"), "literature", Dict(
    "ppp_experiments_commit" => gs.commit,
    "ppp_experiments_branch" => gs.branch,
    "ppp_experiments_dirty"  => gs.dirty,       # scoped to SOURCE_FILES
    "source_files"           => SOURCE_FILES,
    "nadp_seed"              => NADP_SEED,
    "nadp_draws"             => NADP_DRAWS,
    "generated"              => Libc.strftime("%Y-%m-%d", time()),
    "rows"                   => nrow(tbl),
))
println("wrote $(nrow(tbl)) rows -> $csv_out")
println(combine(groupby(tbl, :metabolite), nrow => :n))
```

- [ ] **Step 5: Commit the driver first, then run it (provenance must name a committed source)**

```bash
cd /home/james/projects/PPP_Experiments && git branch --show-current   # feat/phasemap-literature-export
```
Add to `scripts/README.md`, in the grouped-dirs list (alphabetical position), the entry:
```markdown
- `phasemap_export/` — exports `GlycolysisPentosePathwayPhaseMap.jl`'s
  `data/literature_pools.csv` (metabolite-panel literature on the model basis) and its
  `[literature]` provenance: `julia --project scripts/phasemap_export/export_literature_pools.jl`.
  Tests: `julia --project scripts/phasemap_export/test_literature_export.jl`.
```
(Match the surrounding list's exact bullet style when inserting.)
```bash
git add scripts/phasemap_export/export_literature_pools.jl scripts/README.md
git commit -F - <<'EOF'
feat(phasemap-export): add the driver that refreshes the phase-map literature table

Adds the script that writes the phase-map package's literature table and
records which commit of this repo produced it.

- export_literature_pools.jl: import_data.jl -> literature_pools_table over
  the package's panel order -> <phasemap>/data/literature_pools.csv, plus the
  [literature] provenance table (commit, branch, source-scoped dirty flag,
  NADP seed/draws).
- scripts/README.md entry.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KBx8vahm6QPQeQjsQuYbtw
EOF
julia --project scripts/phasemap_export/export_literature_pools.jl
```
Expected: `wrote N rows -> .../GlycolysisPentosePathwayPhaseMap.jl/data/literature_pools.csv` and a per-metabolite count table with 15 metabolites (PGLn absent), NADP = 1000, NADPH = 2. If it throws `ArgumentError` on a value or on non-ASCII text, STOP and report the offending rows (a data issue in `data/Metabolite_Concentrations.csv`, not something to patch in the exporter).
Then check: `grep -A9 '^\[literature\]' /home/james/projects/GlycolysisPentosePathwayPhaseMap.jl/data/provenance.toml` shows `ppp_experiments_dirty = false` and the commit just made.

- [ ] **Step 6: Delete the old table and run the PhaseMap suite**

```bash
cd /home/james/projects/GlycolysisPentosePathwayPhaseMap.jl && git rm -q data/Metabolite_Concentrations.csv
julia --project=. -e 'using Pkg; Pkg.test()'
```
Expected: all testsets PASS (including the shipped-table tests from Step 1).

- [ ] **Step 7: Commit (PhaseMap)**

```bash
cd /home/james/projects/GlycolysisPentosePathwayPhaseMap.jl && git branch --show-current   # feat/ppp-main-resync
git add src/metab_panel.jl src/GlycolysisPentosePathwayPhaseMap.jl test/test_metab_panel.jl test/test_explorer.jl data/literature_pools.csv data/provenance.toml
git commit -F - <<'EOF'
feat(metab-panel): compare model pools against literature on the model's basis

The metabolite panel's grey literature boxes were plotted as published
(mostly whole-cell totals, total NADP(H)) next to model points that are
cytosolic-water, free-NADP(H) concentrations. The boxes now come from a
table exported by PPP_Experiments that is already on the model's basis, so
the two are directly comparable. Free NADP+ has no direct measurement; its
box is derived as free NADPH divided by the measured NADPH/NADP+ ratio.

- data/literature_pools.csv (long format) replaces Metabolite_Concentrations.csv;
  [literature] provenance in data/provenance.toml.
- load_metab_lit reads the long format; mirrored GAP/DHAP corrections removed
  (the exporter applies them).
- y-label states the basis (_METAB_YLABEL).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KBx8vahm6QPQeQjsQuYbtw
EOF
```

---

### Task 5: Regenerate the shipped grids (PhaseMap)

**Files:**
- Data: `data/atp_demand_grid.csv`, `data/demand_grid.csv`, `data/provenance.toml` (`[grids]`)
- Test: `test/test_grid_io.jl`

**Interfaces:**
- Consumes: `regenerate_grid` (Task 2), `G.read_provenance`, `G.load_grid`, `G._PPP_FLUX_KEYS`, `G._PANEL_METAB_ORDER`.
- Produces: shipped grids + `[grids]` table that Task 6's docs and screenshot read.

- [ ] **Step 1: Gate on the PPP.jl checkout**

Run the three commands from Task 2 Step 1. Expected: `fcccae8a196ce9cfba23fae91ca77c4be4540807`, empty, `main`. Otherwise STOP.

- [ ] **Step 2: Write the failing shipped-data tests**

Append to `test/test_grid_io.jl` (and add `CSV` to its `using` line: `using Test, CSV, DataFrames, GlycolysisPentosePathwayPhaseMap`):
```julia
@testset "shipped grids carry the viewer's columns" begin
    for path in (G.default_atp_grid(), G.default_demand_grid())
        cols = names(CSV.read(path, DataFrame; limit = 1))
        for k in G._PPP_FLUX_KEYS
            @test "$(k)_uM_min" in cols
        end
        for s in G._PANEL_METAB_ORDER
            @test "$(s)_uM" in cols
        end
        @test "retcode" in cols && "cycle_index" in cols && "mode" in cols
    end
end

@testset "shipped provenance matches shipped data" begin
    p = G.read_provenance()
    @test p["grids"]["ppp_commit"] == "fcccae8a196ce9cfba23fae91ca77c4be4540807"
    @test p["grids"]["ppp_dirty"] == false
    @test p["grids"]["atp_grid_rows"] == nrow(G.load_grid(G.default_atp_grid()))
    @test p["grids"]["demand_grid_rows"] == nrow(G.load_grid(G.default_demand_grid()))
    @test p["grids"]["atp_grid_nonconverged"] ==
          count(!=(:Terminated), G.load_grid(G.default_atp_grid()).retcode)
    @test !isempty(p["literature"]["ppp_experiments_commit"])
end
```

- [ ] **Step 3: Run to verify failure**

Run: `cd /home/james/projects/GlycolysisPentosePathwayPhaseMap.jl && julia --project=. -e 'using GlycolysisPentosePathwayPhaseMap, Test; include("test/test_grid_io.jl")'`
Expected: the provenance testset ERRORs with `KeyError: key "grids" not found`.

- [ ] **Step 4: Regenerate (background)**

Run with `run_in_background: true` and wait for the completion notification:
```bash
cd /home/james/projects/GlycolysisPentosePathwayPhaseMap.jl && julia --project=/home/james/projects/PPP_Experiments -e 'using GlycolysisPentosePathwayPhaseMap, PentosePhosphatePathway; regenerate_grid()' 2>&1 | tail -20
```
Expected: `Regenerated grids … atp_rows = 1584 demand_rows = 144`. Then:
```bash
grep -A20 '^\[grids\]' data/provenance.toml
```
Expected: `ppp_commit = "fcccae8a196ce9cfba23fae91ca77c4be4540807"`, `ppp_branch = "main"`, `ppp_dirty = false`, and the `[literature]` table still present.

- [ ] **Step 5: Run the full suite**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS. (The explorer test needs ≥1 converged cell and 11 ATP levels; both hold if the run completed.)

- [ ] **Step 6: Commit**

```bash
git branch --show-current   # feat/ppp-main-resync
git add data/atp_demand_grid.csv data/demand_grid.csv data/provenance.toml test/test_grid_io.jl
git commit -F - <<'EOF'
data(grids): regenerate both demand grids from PentosePhosphatePathway.jl main

Recomputes every precomputed operating point the explorer and phase map show,
using the current kinetic model (47 commits newer than the previous grids:
new G6PD rate law, corrected starting concentrations, free NADP(H)/NAD(H)
pools, an ATP-consuming R5P demand sink, and other parameter updates).

- data/atp_demand_grid.csv (1584 rows) and data/demand_grid.csv (144 rows)
  from PPP.jl main @ fcccae8 (clean), same axes and solver settings.
- [grids] provenance table; tests pin provenance to the shipped files.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KBx8vahm6QPQeQjsQuYbtw
EOF
```

---

### Task 6: Docs, screenshot, version, local CLAUDE.md (PhaseMap)

**Files:**
- Modify: `docs/src/recompute.md`, `docs/src/tutorial.md`, `docs/src/index.md`, `Project.toml`, `docs/src/assets/explorer.png`, `README.md` (only if Step 1 finds a contradicted claim)
- Local only: `CLAUDE.md`

**Interfaces:**
- Consumes: regenerated grids (Task 5), `literature_pools.csv` (Task 4).

- [ ] **Step 1: Check result claims against the new grids**

Run:
```bash
cd /home/james/projects/GlycolysisPentosePathwayPhaseMap.jl && julia --project=. -e '
using GlycolysisPentosePathwayPhaseMap, DataFrames
const G = GlycolysisPentosePathwayPhaseMap
for p in (G.default_demand_grid(), G.default_atp_grid())
    df = G.load_grid(p)
    println(basename(p), ": ", sort(collect(pairs(Dict(m => count(==(m), df.mode) for m in unique(df.mode))))),
            "  nonconverged=", count(!=(:Terminated), df.retcode))
end'
grep -n -i "blue\|red\|gray\|reverse\|linear\|cycle corner\|sign-flip\|contour" README.md docs/src/tutorial.md docs/src/science.md docs/src/index.md
```
For each grep hit that states what the grids *show* (e.g. "red cells = reverse PPP", "a black contour marks the sign-flip boundary", `phasemap` `reps` defaulting to a cycle and a reverse corner): if the mode counts show that region is absent (count 0) in the 2-D grid, edit that sentence to describe what the grid now shows; otherwise leave it. Definitions (what the index means, what a mode is) are not result claims and stay. Also run `julia --project=. -e 'using GlycolysisPentosePathwayPhaseMap; phasemap(save=tempname()*".png")'` — Expected: no error (it needs one cycle and one reverse representative cell; if it errors because a mode is absent, STOP and report rather than changing `phasemap`'s defaults).

- [ ] **Step 2: Update `docs/src/index.md`**

Replace "is a visualization package for the red-blood-cell\nglycolysis + pentose-phosphate-pathway (PPP)" with "is a visualization package for the general mammalian\nglycolysis + pentose-phosphate-pathway (PPP)" (keep the line wrap ≤ 100 chars).

- [ ] **Step 3: Update the tutorial's metabolite-panel paragraph**

In `docs/src/tutorial.md`, replace the paragraph starting "Below the network sits a full-width **metabolite box-plot panel**" through "flag known model discrepancies." with:
```markdown
Below the network sits a full-width **metabolite box-plot panel** covering 16 metabolites
(`Glucose, G6P, F6P, F16BP, DHAP, GAP, ATP, NADP, NADPH, PGLn, PGA, Ru5P, R5P, X5P, S7P, E4P`).
For the currently selected cell, each metabolite's steady-state pool (µM, log y-axis) is drawn as
a **purple diamond** over a **static literature box-and-whisker** from cultured mammalian cells.
Both are on the model's pool basis: µM of cytosolic water, with NADP⁺ and NADPH as **free**
(unbound) concentrations and GAP/DHAP as the free aldehyde/ketone. The literature values were
converted to that basis before shipping (see [Recompute](@ref)); the boxes include drug-treated
and glucose-pulse measurements as well as basal ones. Free NADP⁺ has no direct measurement, so
its box is **derived** — free NADPH divided by the measured free NADPH/NADP⁺ ratio — and PGLn
has no measurements, so it has no box. The diamonds redraw on every click and slider move; the
boxes are static. Diamonds sitting off their boxes flag known model discrepancies.
```

- [ ] **Step 4: Rewrite `docs/src/recompute.md`**

Replace the whole file with:
````markdown
# Recompute

The package ships **precomputed** data in its `data/` directory, so `explore()` and
`phasemap()` work instantly with no kinetic model and no ODE solve:

- `atp_demand_grid.csv`, `demand_grid.csv` — steady-state fluxes and pools over the demand grids;
- `literature_pools.csv` — the metabolite panel's literature values, on the model's pool basis;
- `provenance.toml` — which builds produced them: `[grids]` and `[literature]`.

Each data source is refreshed separately, and each refresh rewrites only its own provenance table.

## Regenerating the grids

Regenerating requires the kinetic model and solves the glycolysis + PPP ODE system at every
demand cell. It is gated behind a package extension, so the model and solver stack are **not**
dependencies of the visualization package.

1. Install the kinetic model alongside this package:

   ```julia
   using Pkg
   Pkg.add(url="https://github.com/DenisTitovLab/PentosePhosphatePathway.jl")
   ```

2. Load it. This activates the recompute extension (`PentosePhosphatePathwayExt`):

   ```julia
   using GlycolysisPentosePathwayPhaseMap
   using PentosePhosphatePathway
   ```

3. Regenerate:

   ```julia
   regenerate_grid()   # writes data/*.csv and the [grids] table of data/provenance.toml
   ```

`[grids]` records the PentosePhosphatePathway.jl commit (with branch and a dirty flag when it is a
git checkout; otherwise its version), the CellMetabolism version, the grid axes, and how many
cells did not reach steady state (those render gray).

If you call `regenerate_grid` **without** the extension loaded, it raises an error telling you
to `using PentosePhosphatePathway` first.

!!! note
    `regenerate_grid` solves the kinetic ODE model at every demand cell and takes **minutes**.
    The first `using` / solve in a session also pays Julia's compilation cost. This is expected,
    not a hang.

## Refreshing the literature table

`literature_pools.csv` is exported from the
[PPP_Experiments](https://github.com/DenisTitovLab/PPP_Experiments) data repository, where the
measurements and their basis conversion live (each row's reported basis → µM of cytosolic water
by species distribution class; free-pool NADP(H); GAP/DHAP free-form corrections). From a
PPP_Experiments checkout next to this package:

```sh
julia --project scripts/phasemap_export/export_literature_pools.jl
```

This rewrites `data/literature_pools.csv` and the `[literature]` table (PPP_Experiments commit,
source files, and the seed and draw count of the derived free-NADP⁺ values).
````

- [ ] **Step 5: Re-render the screenshot and bump the version**

```bash
cd /home/james/projects/GlycolysisPentosePathwayPhaseMap.jl && julia --project=. scripts/render_readme_screenshot.jl
```
Expected: `wrote …/docs/src/assets/explorer.png`. Open it with the Read tool and confirm the panel shows boxes, diamonds, and the new y-label.
In `Project.toml` change `version = "0.1.0"` to `version = "0.2.0"`.

- [ ] **Step 6: Build the docs**

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```
Expected: `makedocs` finishes without errors (`deploydocs` skips outside CI). `docs/build/` is gitignored.

- [ ] **Step 7: Rewrite the local CLAUDE.md (not committed)**

Overwrite `/home/james/projects/GlycolysisPentosePathwayPhaseMap.jl/CLAUDE.md` with:
````markdown
# CLAUDE.md

Guidance for Claude Code in this repo. This file is gitignored (local only).

## What this package is

Visualization of the **general mammalian** glycolysis + pentose-phosphate-pathway (PPP) *demand
phase map*, built on the kinetic model `PentosePhosphatePathway.jl` (PPP.jl). It ships
precomputed data in `data/`, so the viewers need no kinetic model and no ODE solve:

- `atp_demand_grid.csv` (3-D ATP × NADPH × R5P) and `demand_grid.csv` (2-D NADPH × R5P):
  steady-state fluxes (`*_uM_min`) and pools (`*_uM`) per demand cell.
- `literature_pools.csv`: metabolite-panel literature values, already on the model's basis.
- `provenance.toml`: `[grids]` (PPP.jl commit/branch/dirty, axes, non-converged counts) and
  `[literature]` (PPP_Experiments commit). Each writer replaces only its own table.

Public API (`src/GlycolysisPentosePathwayPhaseMap.jl`): `explore()` (GLMakie window),
`phasemap(; save)` (CairoMakie figure), `regenerate_grid(; …)` (needs the extension).

## Commands

```sh
julia --project=. -e 'using Pkg; Pkg.test()'                       # core suite (no PPP.jl)
julia --project=. -e 'using GlycolysisPentosePathwayPhaseMap, Test; include("test/test_grid_io.jl")'
julia --project=/home/james/projects/PPP_Experiments test/manual/test_extension.jl   # extension
julia --project=docs docs/make.jl                                   # docs
julia --project=. scripts/render_readme_screenshot.jl               # docs/src/assets/explorer.png
```
GLMakie needs an X display even to load; on a headless box prefix with `xvfb-run -a`.

## Refreshing shipped data (two independent steps)

1. **Grids** — from the PPP_Experiments env (its `[sources]` `develop` this package and
   `../PentosePhosphatePathway.jl` by path). First confirm the PPP.jl checkout is `main` and clean:
   `julia --project=/home/james/projects/PPP_Experiments -e 'using GlycolysisPentosePathwayPhaseMap, PentosePhosphatePathway; regenerate_grid()'`
   (minutes; run in the background).
2. **Literature** — from the PPP_Experiments root:
   `julia --project scripts/phasemap_export/export_literature_pools.jl`.
   All basis conversion (row `Basis` × species class → cytosolic water; free NADP(H); GAP/DHAP
   free-form) lives there. Never re-implement it here.

Commit the regenerated CSVs together with `data/provenance.toml`.

## Architecture

- `src/grid_io.jl` — CSV contract shared by renderers and the extension (`load_grid`,
  `row_to_fluxes`, `row_to_pools`, `pentose_cycle_index`, `classify_mode`, `_PANEL_METAB_ORDER`).
- `src/provenance.jl` — `read_provenance`, `update_provenance!`, `git_state`.
- `src/flux_map.jl` — network schematic; `GAP`/`F6P` are drawn twice as a display device only.
- `src/metab_panel.jl` — literature boxes + model diamonds; `load_metab_lit` reads the long table.
- `src/phasemap.jl`, `src/explorer.jl` — static figure and interactive window
  (`build_explorer(df)` is display-free and offscreen-testable).
- `ext/PentosePhosphatePathwayExt.jl` — the ODE sweep (warm-started down the NADPH axis).

**Extension hook pattern (keep it):** an extension may not define a method on a parent-owned
function (precompile error). The core declares `const _REGEN_HOOK = Ref{Any}(nothing)`; the
extension's `__init__` sets it; `regenerate_grid` calls through it. The extension triggers on all
six weakdeps; `using PentosePhosphatePathway` loads the other five.

**Backends:** `import GLMakie, CairoMakie` (never `using` both: re-export clash). `explore`
activates GLMakie, `phasemap` CairoMakie.

## Model conventions

- Demand levers as fractions φ of supply ceilings anchored on `RPI_cap = RPI_Vmax·RPI_Conc`:
  NADPH `6·RPI_cap`, R5P `1.5·RPI_cap`, ATP `2·HK1_cap`. Grids are log-spaced in φ.
- Pentose-cycle index `(V_ox − V_R5Pase)/V_ox`: ≈1 cycle, ≈0 linear, <0 reverse.
- Initial state is `copy(PPP.glycolysis_ppp_init_conc)` — never hand-merge PPP and CellMetabolism
  inits (that skips `PPP_glycolysis_init_overrides`). Params: `merge_LArrays(PPP_params,
  CellMetabolism.glycolysis_params)` (PPP.jl wins).
- Solver: `Rodas5P()`, `TerminateSteadyState(1e-13, 1e-8)`, `abstol = 1e-14`, `reltol = 1e-9`.
- Non-converged cells keep `retcode != :Terminated` and render as the colormap's `nan_color`.
- 6-phosphogluconate is `PGA`, never `6PG`.

## Conventions

- Docstrings bind to the next expression: no `const` between a docstring and its function
  (`checkdocs = :exports` fails otherwise).
- CI runs tests under `xvfb-run`; the docs job starts a persistent Xvfb.
- Section banners: plain ASCII `#`, 90 characters wide, title centred.
- Registration in General is blocked until PPP.jl is registered (it is a weakdep).
````
Then confirm: `wc -l CLAUDE.md` → < 100, and `git status --porcelain CLAUDE.md` → empty (ignored).

- [ ] **Step 8: Commit**

```bash
git branch --show-current   # feat/ppp-main-resync
git add docs/src/recompute.md docs/src/tutorial.md docs/src/index.md docs/src/assets/explorer.png Project.toml
# plus README.md / docs/src/science.md only if Step 1 edited them
git commit -F - <<'EOF'
docs: describe the two-step data refresh and the panel's literature basis; bump to 0.2.0

Updates the documentation for the resynced data: how to refresh the grids and
the literature table separately, what the provenance file records, and what
the metabolite panel's literature boxes now contain. Also corrects the home
page's red-blood-cell framing to general mammalian cells and re-renders the
explorer screenshot from the new grids.

- recompute.md: grids via the extension, literature via PPP_Experiments,
  data/provenance.toml.
- tutorial.md: model basis, free NADP(H), derived NADP+, treated conditions,
  no PGLn box.
- explorer.png re-rendered; version 0.1.0 -> 0.2.0 (shipped results change).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01KBx8vahm6QPQeQjsQuYbtw
EOF
```

---

### Task 7: Final verification and hand-off

- [ ] **Step 1: Full verification**

```bash
cd /home/james/projects/GlycolysisPentosePathwayPhaseMap.jl && julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=/home/james/projects/PPP_Experiments test/manual/test_extension.jl
cd /home/james/projects/PPP_Experiments && julia --project scripts/phasemap_export/test_literature_export.jl
git -C /home/james/projects/PentosePhosphatePathway.jl branch --show-current   # main
```
Expected: all PASS; PPP.jl on `main`. Note: the manual extension test writes only to a temp dir.

- [ ] **Step 2: Integrate**

Use superpowers:finishing-a-development-branch for both branches. PhaseMap `feat/ppp-main-resync` → push and open a PR to `main` (already approved). Push with the sandbox disabled if `.git/config.lock` blocks `push -u`. For PPP_Experiments `feat/phasemap-literature-export`, ask James whether to open a PR or merge locally.
