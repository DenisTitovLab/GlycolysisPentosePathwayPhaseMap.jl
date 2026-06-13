# GlycolysisPentosePathwayPhaseMap.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://DenisTitovLab.github.io/GlycolysisPentosePathwayPhaseMap.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://DenisTitovLab.github.io/GlycolysisPentosePathwayPhaseMap.jl/dev/)
[![Build Status](https://github.com/DenisTitovLab/GlycolysisPentosePathwayPhaseMap.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DenisTitovLab/GlycolysisPentosePathwayPhaseMap.jl/actions/workflows/CI.yml?query=branch%3Amain)

A visualization package for the mammalian glycolysis + pentose-phosphate-pathway (PPP)
*demand phase map*. It shows how perturbing metabolic demand — for NADPH, ribose-5-phosphate
(R5P), and ATP — reorganizes the **direction and scale** of PPP flux, switching the pathway
between linear, full pentose-cycle, and reverse (non-oxidative ribose-producing) operation. The
package ships precomputed grids of steady-state fluxes and metabolite pools, so its tools work
**instantly, with no kinetic model and no ODE solve**.

![The explore() interactive demand phase map: pentose-cycle-index heatmap, glycolysis + PPP flux
network, and steady-state metabolite panel for the selected cell.](docs/src/assets/explorer.png)

The `explore()` window: click any cell on the pentose-cycle-index heatmap (left) to redraw the
glycolysis + PPP flux network (right) and the steady-state metabolite pools (bottom); the slider
re-slices the cached 3-D grid by ATP demand.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/DenisTitovLab/GlycolysisPentosePathwayPhaseMap.jl")
```

## Quickstart

```julia
using GlycolysisPentosePathwayPhaseMap
explore()                        # interactive GLMakie window (cached grid, no ODE solve)
phasemap(save="phasemap.png")    # static Layout-A phase map figure (PNG)
```

## What you're looking at

Three demand levers drive the integrated glycolysis + PPP model, each set as a fraction **φ** of
that species' *sustainable supply ceiling*, so all axes mean "percent of what the pathway can
make." The ceilings are anchored on ribose-5-phosphate isomerase (RPI), the carbon bottleneck:

| Species | Supply ceiling |
|---|---|
| NADPH | `6 · RPI_cap` |
| R5P   | `1.5 · RPI_cap` |
| ATP   | `2 · HK1_cap` |

Each operating point is summarized by the **pentose-cycle index** `(V_ox − V_R5Pase) / V_ox`:
≈ 1 is a full forward pentose cycle, ≈ 0 is linear PPP, and < 0 is reverse (non-oxidative R5P
production). The explorer's heatmap colors cells by this index (blue = cycle, red = reverse,
gray = non-converged) and redraws the glycolysis + PPP flux network for any cell you click; an
ATP-demand slider re-slices the cached 3-D grid live.

## Recompute

The shipped grids are precomputed. To regenerate them you install the kinetic model
(`Pkg.add(url="https://github.com/DenisTitovLab/PentosePhosphatePathway.jl")`), `using
PentosePhosphatePathway` to activate the recompute extension, then call `regenerate_grid(; ...)`
(solves the ODE model; takes minutes). See the
[Recompute docs](https://DenisTitovLab.github.io/GlycolysisPentosePathwayPhaseMap.jl/stable/recompute/).

## Documentation

Full documentation: <https://DenisTitovLab.github.io/GlycolysisPentosePathwayPhaseMap.jl/stable/>
