# Respiration what-if grids (branch `feat/atp-synthase-respiration`)

The phase map re-solved on the model with respiration added as an ATP-producing pseudo-enzyme,
for a 2x2 comparison of respiration x F26BP. Each folder holds an `atp_demand_grid.csv` (3-D),
a `demand_grid.csv` (2-D, ATP demand 0.25) and a `provenance.toml` (`[grids]` only).

| Respiration fraction r | F26BP 5 µM | F26BP 1 µM |
|---|---|---|
| 0 (no respiration) | `data/` (the shipped grids) | `r0_f26bp1uM/` |
| 0.5 | `r0.5_f26bp5uM/` | `r0.5_f26bp1uM/` |

```julia
using GlycolysisPentosePathwayPhaseMap
const G = GlycolysisPentosePathwayPhaseMap
explore(grid = joinpath(G.DATA_DIR, "respiration", "r0.5_f26bp5uM", "atp_demand_grid.csv"))
phasemap(grid = joinpath(G.DATA_DIR, "respiration", "r0.5_f26bp5uM", "demand_grid.csv"))
```

## The model

Construction of Kober et al. 2026 (bioRxiv 10.64898/2026.01.31.703051, Fig 1B/S2; code in
CellMetabolism.jl `examples/2026_F26BP_Figures/`), implemented in `ext/respiration.jl`:

- **Reaction** `ATPSynthase: ADP + Pi -> ATP`, reversible Michaelis–Menten with competitive
  product inhibition. Km_ATP = 1000 M, Km_ADP = Km_Pi = 10 µM, Keq = 1e6 M⁻¹. These constants
  are the paper's model choices, **not measurements**. The reaction has no carbon or NADH cost.
- **Fixed total ATP-supply capacity** C = 0.3 · 2·HK1_cap (the most glycolysis delivers).
  `ATPSynthase_Vmax = r · C`. The 13 glycolytic enzyme concentrations (GLUT1, HK1, GPI, PFKP,
  ALDO, TPI, GAPDH, PGK, PGAM, ENO, PKM2, LDH, MCT) are multiplied by (1 − r). PPP enzymes,
  sinks, kinetic constants and initial conditions are those of PPP.jl `main`.
- **Axes are unchanged**: the ATPase, NADPH and R5P demand ceilings come from the unscaled
  parameters, so each cell is the same absolute demand in all four corners.
- Extra columns: `V_PKM2_uM_min` (glycolytic ATP production), `V_ATPSynthase_uM_min`,
  `V_ATPase_uM_min`, `Phosphate_uM`.

## Regenerate

```sh
julia --project=/home/james/projects/PPP_Experiments scripts/regenerate_respiration_grids.jl
```

Solved 2026-10-15 against PentosePhosphatePathway.jl `main` @ `fcccae8` (a few minutes for all
four solves). The script first re-solves the r = 0, 5 µM corner through the respiration code
path and refuses to write unless it reproduces `data/`. On this run it was bit-identical
(max rel. diff 0) for both grids.

## Results (2026-10-15)

Physiological cell (NADPH φ 0.059, R5P φ 0.1068, ATP demand 0.25 · 2·HK1_cap), from the 2-D
grids, plus non-converged and undetermined (ATP-collapsed) cell counts over the 1584-cell 3-D
grids:

| Corner | Glycolytic ATP (V_PKM2) | Respiration | Resp. share | V_HK1 | Cycle index | ATP | ADP | Non-conv. | Undetermined |
|---|---|---|---|---|---|---|---|---|---|
| r 0, 5 µM | 7376 | 0 | 0% | 3717 | 0.063 | 7472 | 28.8 | 0 | 184 |
| r 0.5, 5 µM | 3942 | 3431 | 46.5% | 2000 | 0.066 | 7462 | 39.1 | 0 | 0 |
| r 0.5, 1 µM | 3621 | 3750 | 50.9% | 1839 | 0.076 | 7437 | 63.8 | 0 | 0 |
| r 0, 1 µM | 7370 | 0 | 0% | 3714 | 0.080 | 7428 | 73.2 | 1 | 193 |

Fluxes are in µM/min and pools in µM, on the model's cytosolic-water basis.

- **Kober's competition is reproduced.** Respiration takes about half of the demand-set ATP
  flux, and glycolysis (and glucose uptake through HK1) falls one-for-one. Lowering F26BP from
  5 to 1 µM moves 4.4 percentage points of ATP supply from glycolysis to respiration, the
  direction Kober reports.
- **The physiological operating point barely moves.** With respiration, the cycle index there
  changes by 0.003 (5 µM) and 0.004 (1 µM), even though glucose uptake roughly halves.
- **Respiration removes the high-NADPH collapse.** With r = 0, 16–17 cells per ATP level
  collapse (ATP → 0, R5P sink running backward), at every ATP demand: the whole top NADPH row
  (φ = 1) plus the high-R5P corner of the next two rows (the gray region of the static phase
  map). With r = 0.5, all 1584 cells keep a defined cycle index, and cells at NADPH φ ≥ 0.29
  hold a median ATP of 7.4 mM, against 0.86 mM with r = 0. This rescue relies on respiration
  needing no glycolytic carbon. In a cell, mitochondria burn pyruvate, which a collapsed
  glycolysis would no longer supply. Read it as "an ATP source independent of glycolysis
  prevents the collapse", not as a prediction for real respiring cells.
