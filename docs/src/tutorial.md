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
- **gray** cells = operating points with no meaningful index: non-converged, oxidative flux ≈ 0,
  or the R5P sink running backward after ATP has collapsed (masked, never silently dropped).

A black contour marks the sign-flip boundary between the cycle and reverse regions.

### Clicking cells

Click any cell in the heatmap. The large combined glycolysis + PPP **flux network** on the right
redraws from that operating point's cached fluxes. The footer reports the absolute `max |flux|`
(µM/min) for that cell, its mode, and its cycle index.

### The ATP-demand slider

Below the panels is an **ATP-demand slider** that varies ATP demand across log-spaced levels
(1 %–25 % of maximal ATP supply, i.e. 0.02–0.5 × HK1 capacity). Moving it re-slices the cached
3-D grid: the heatmap, its sign-flip contour, and the flux network all update instantly — again,
with no ODE re-solve. The explorer opens at the 25 % level (ATP demand 0.5 × HK1 capacity) on the
cell nearest **physiological demand**: PentosePhosphatePathway.jl's default NADPH and R5P sinks,
with the glycolytic activator F26BP at 5 µM (PentosePhosphatePathway.jl's own default is 0).

### The flux network

The network is the glycolysis spine plus the oxidative branch (G6PD → PGLS → PGD) and the
reversible non-oxidative PPP (RPI, RPE, TKT, TALDO1). Each arrow's **color encodes direction**
(forward vs reverse) and its **width is proportional to |flux|**, autoscaled per selected cell.
Metabolite nodes carry their model symbols (`G6P`, `F16BP`, `Ru5P`, …); arrows carry their
enzyme names.

### The metabolite panel

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

## The static phase map

For a publication-style figure, render the static Layout-A phase map with CairoMakie:

```julia
phasemap()                       # returns a Figure
phasemap(save="phasemap.png")    # returns the Figure and writes a PNG
```

The figure is a cycle-index heatmap over the NADPH × R5P demand plane (with the sign-flip
contour) plus two representative glycolysis + PPP flux schematics — the cell nearest
physiological NADPH and R5P demand, and the most strongly reverse cell. The 2-D grid is solved at
an ATP demand of 25 % of maximal supply (0.5 × HK1 capacity).
