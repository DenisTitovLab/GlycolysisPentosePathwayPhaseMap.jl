```@meta
CurrentModule = GlycolysisPentosePathwayPhaseMap
```

# GlycolysisPentosePathwayPhaseMap.jl

`GlycolysisPentosePathwayPhaseMap.jl` is a visualization package for the red-blood-cell
glycolysis + pentose-phosphate-pathway (PPP) *demand phase map*. It shows how perturbing
metabolic demand — for NADPH, for ribose-5-phosphate (R5P), and for ATP — reorganizes the
**direction and scale** of flux through the PPP, switching the pathway between linear, full
pentose-cycle, and reverse (non-oxidative ribose-producing) operation. The package ships
precomputed grids of steady-state fluxes and metabolite pools, so the interactive explorer and
the static phase map work **instantly, with no kinetic model and no ODE solve**. It recapitulates
the non-oxidative-PPP / NADPH-homeostasis phenomenology of Feng et al., *PNAS* 2026,
123(8):e2526325123.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/DenisTitovLab/GlycolysisPentosePathwayPhaseMap.jl")
```

## Quickstart

Open the interactive explorer (a native GLMakie window):

```julia
using GlycolysisPentosePathwayPhaseMap
explore()
```

Or render the static Layout-A phase map figure:

```julia
phasemap()                       # returns a Figure
phasemap(save="phasemap.png")    # also writes a PNG
```

Both read precomputed data shipped in the package's `data/` directory, so they return
immediately without solving the kinetic model.

See the [Tutorial](@ref) for a guided walk-through, [The science](@ref) for the
supply/demand convention and cycle-index definition, and [Recompute](@ref) for regenerating the
cached grids from the kinetic model.
