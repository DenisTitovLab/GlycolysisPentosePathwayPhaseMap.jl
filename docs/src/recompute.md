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
[PPP_Experiments](https://github.com/JamesMbata/PPP_Experiments) data repository, where the
measurements and their basis conversion live (each row's reported basis → µM of cytosolic water
by species distribution class; free-pool NADP(H); GAP/DHAP free-form corrections). From a
PPP_Experiments checkout next to this package:

```sh
julia --project scripts/phasemap_export/export_literature_pools.jl
```

This rewrites `data/literature_pools.csv` and the `[literature]` table (PPP_Experiments commit,
source files, and the seed and draw count of the derived free-NADP⁺ values).
