# Tutorial

This page walks through the interactive explorer and the static phase map. Everything below
reads precomputed data shipped with the package — no kinetic model, no ODE solve.

## Launching the explorer

```julia
using GlycolysisPentosePathwayPhaseMap
explore()
```

`explore()` opens a native GLMakie window over a cached 3-D ATP × NADPH × R5P demand grid. It
blocks until you close the window. To use a different grid CSV, pass a path:
`explore("path/to/grid.csv")`.

<!-- A screenshot of the explorer will be added here (assets/explorer.png). -->

### Reading the cycle-index heatmap

The left panel is a heatmap of the **pentose-cycle index** over the NADPH × R5P demand plane
(see [The science](@ref) for the definition). Under the diverging colormap:

- **blue** cells = forward pentose cycle (index ≈ 1),
- **red** cells = reverse PPP (index < 0, non-oxidative R5P production),
- **gray** cells = non-converged operating points (masked, never silently dropped).

A black contour marks the sign-flip boundary between the cycle and reverse regions.

### Clicking cells

Click any cell in the heatmap. The large combined glycolysis + PPP **flux network** on the right
redraws from that operating point's cached fluxes. The footer reports the absolute `max |flux|`
(µM/min) for that cell, its mode, and its cycle index.

### The ATP-demand slider

Below the panels is an **ATP-demand slider** that varies ATP demand across log-spaced levels
(roughly 1 %–20 % of maximal ATP supply). Moving it re-slices the cached 3-D grid: the heatmap,
its sign-flip contour, and the flux network all update instantly — again, with no ODE re-solve.
The explorer opens at the level nearest 10 % of ATP supply.

### The flux network

The network is the glycolysis spine plus the oxidative branch (G6PD → PGLS → PGD) and the
reversible non-oxidative PPP (RPI, RPE, TKT, TALDO1). Each arrow's **color encodes direction**
(forward vs reverse) and its **width is proportional to |flux|**, autoscaled per selected cell.
Metabolite nodes carry their model symbols (`G6P`, `F16BP`, `Ru5P`, …); arrows carry their
enzyme names.

### The metabolite panel

Below the network sits a full-width **metabolite box-plot panel** covering 16 metabolites
(`Glucose, G6P, F6P, F16BP, DHAP, GAP, ATP, NADP, NADPH, PGLn, PGA, Ru5P, R5P, X5P, S7P, E4P`).
For the currently selected cell, each metabolite's steady-state pool (µM, log y-axis) is drawn as
a **purple diamond** over a **static literature box-and-whisker** from measured RBC concentrations.
The diamonds redraw on every click and slider move; the boxes are static. Diamonds sitting off
their boxes flag known model discrepancies.

## The static phase map

For a publication-style figure, render the static Layout-A phase map with CairoMakie:

```julia
phasemap()                       # returns a Figure
phasemap(save="phasemap.png")    # returns the Figure and writes a PNG
```

The figure is a cycle-index heatmap over the NADPH × R5P demand plane (with the sign-flip
contour) plus two representative glycolysis + PPP flux schematics — one for a cycle operating
point and one for a reverse operating point.
