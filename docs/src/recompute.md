# Recompute

The package ships **precomputed** grid CSVs in its `data/` directory. That is why `explore()` and
`phasemap()` work instantly with no kinetic model installed and no ODE solve — they only read and
visualize cached steady-state fluxes and metabolite pools.

Regenerating those grids is a separate, heavier operation: it requires the kinetic model and
solves the glycolysis + PPP ODE system at every demand cell. This is gated behind a package
extension so the model and solver stack are **not** dependencies of the visualization package.

## How to recompute

1. Install the kinetic model alongside this package:

   ```julia
   using Pkg
   Pkg.add(url="https://github.com/DenisTitovLab/PentosePhosphatePathway.jl")
   ```

2. Load it. This activates the recompute extension (`PentosePhosphatePathwayExt`), which pulls in
   the solver stack and wires up `regenerate_grid`:

   ```julia
   using GlycolysisPentosePathwayPhaseMap
   using PentosePhosphatePathway
   ```

3. Regenerate the grid CSVs:

   ```julia
   regenerate_grid(; ...)
   ```

If you call `regenerate_grid` **without** the extension loaded, it raises a friendly error telling
you to `using PentosePhosphatePathway` first.

!!! note
    `regenerate_grid` solves the kinetic ODE model at every demand cell and takes **minutes**.
    The first `using` / solve in a session also pays Julia's compilation cost. This is expected,
    not a hang.
