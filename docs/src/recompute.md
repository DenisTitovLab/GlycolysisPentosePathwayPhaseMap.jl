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

`regenerate_grid(; f26bp_uM = 5.0)` sets the F26BP value (µM of cytosolic water); the default is the phase map's physiological 5 µM, not PentosePhosphatePathway.jl's default of 0.

`nadph_phis` and `r5p_phis` set the demand axes as strictly increasing fractions of each supply
ceiling, for example `regenerate_grid(; nadph_phis = [0.01, 0.1, 0.5, 1.0])`; the viewer draws one
equal-sized cell per level. The defaults are the axes described in [The science](@ref). Grids
written before 0.4.0 sit on 12 evenly log-spaced levels (φ = 10⁻³ to 1) and record their axes as
`phi_lo`/`phi_hi`/`n_nadph`/`n_r5p` in `[grids]`; they still render.

Grids written before 0.3.0 lack the cofactor pool columns. `explore()` rejects them with an error
asking you to regenerate, while `phasemap()` still renders them.

Every cell is solved from PentosePhosphatePathway.jl's initial state, so its steady state does
not depend on the order the sweep visits cells. Fluxes are cached in µM/min and pools in µM.

`[grids]` records the PentosePhosphatePathway.jl commit (with branch and a dirty flag when it is a
git checkout; otherwise its version), the CellMetabolism version, the grid axes, the F26BP value, and how many
cells did not reach steady state (those render gray).

If you call `regenerate_grid` **without** the extension loaded, it raises an error telling you
to `using PentosePhosphatePathway` first.

!!! note
    `regenerate_grid` solves the kinetic ODE model at every demand cell and takes **minutes**.
    The first `using` / solve in a session also pays Julia's compilation cost. This is expected,
    not a hang.

## Refreshing the literature table

`literature_pools.csv` is exported from the
[PPP_Experiments](https://github.com/JamesMbata/PPP_Experiments) data repository, where the
measurements and their basis conversion live (each row's reported basis → µM of cytosolic water
by species distribution class; free-pool cofactors (ADP, AMP, NAD⁺, NADH, NADP⁺, NADPH: free-pool
measurements and same-cell-line matched pairs only); GAP/DHAP free-form corrections). From a
PPP_Experiments checkout next to this package:

```sh
julia --project scripts/phasemap_export/export_literature_pools.jl
```

This rewrites `data/literature_pools.csv` and the `[literature]` table (PPP_Experiments commit and
source files). Each row's `method` column (`metabolomics`, `sensor` or `equilibrium`) sets its
marker in the panel. The export refuses to run unless `_FREE_POOL_PANEL` equals PPP_Experiments'
`FREE_POOL_METABS` (same order) and every `_FREE_EQ_TOTAL_PANEL` entry is in
`FREE_EQ_TOTAL_METABS`.
