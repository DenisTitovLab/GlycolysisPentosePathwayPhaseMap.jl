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

Implemented in `ext/respiration.jl`. The capacity construction is that of Kober et al. 2026
(bioRxiv 10.64898/2026.01.31.703051, Fig 1B/S2; code in CellMetabolism.jl
`examples/2026_F26BP_Figures/`), but respiration here burns glycolytic carbon:

- **Reaction** `ATPSynthase: Pyruvate + NADH + 15 ADP + 15 Pi -> NAD + 3 CO2 + 15 ATP`.
  - Full pyruvate oxidation gives 12.5 ATP. The cytosolic NADH carried in by the
    malate–aspartate shuttle gives 2.5 more. With glycolysis that totals about 32 ATP per
    glucose.
  - Consuming the NADH is required. The model conserves BPG + 3PG + 2PG + PEP + Pyr − NADH, so
    a pyruvate-only sink could carry no steady-state flux.
  - O2 and H2O are implicit. CO2 is the model's fixed boundary species.
- **Rate law**: irreversible respiratory control, per pyruvate.
  - It is Kober's ADP·Pi saturation term (Km_ADP = Km_Pi = 10 µM, Km_ATP = 1000 M) ×
    Pyr/(140 µM + Pyr) × NADH/(0.02 µM + NADH).
  - The pyruvate constant is LDH's Km.
  - The ADP, Pi and NADH constants are **not measurements**. The NADH constant sits about 10×
    below the model's free NADH (0.22 µM), so respiration is ADP-controlled at normal redox and
    stops only when NADH or pyruvate drain.
- **Fixed total ATP-supply capacity** C = 0.3 · 2·HK1_cap (the most glycolysis delivers).
  - `ATPSynthase_Vmax = r · C / 15` (per pyruvate), so respiratory ATP capacity is r · C.
  - The 13 glycolytic enzyme concentrations (GLUT1, HK1, GPI, PFKP, ALDO, TPI, GAPDH, PGK,
    PGAM, ENO, PKM2, LDH, MCT) are multiplied by (1 − r).
  - PPP enzymes, sinks, kinetic constants and initial conditions are those of PPP.jl `main`.
- **Axes are unchanged**: the ATPase, NADPH and R5P demand ceilings come from the unscaled
  parameters, so each cell is the same absolute demand in all four corners.
- **Extra columns**: `V_PKM2_uM_min` (glycolytic ATP production), `V_ATPSynthase_uM_min` (per
  pyruvate; respiratory ATP is 15× it), `V_LDH_uM_min`, `V_ATPase_uM_min`, `Phosphate_uM`,
  `Pyruvate_uM`. At steady state V_GAPDH = V_LDH + V_ATPSynthase (NADH balance).

The first version of this branch (`24cb571`) used Kober's carbon-free `ADP + Pi -> ATP`. Its
results are kept below for comparison.

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
grids. Fluxes are in µM/min, with respiration in ATP units. Pools are in µM on the model's
cytosolic-water basis. "Non-conv." counts the 1584-cell 3-D grid.

| Corner | Glycolytic ATP (V_PKM2) | Respiratory ATP | Resp. share | V_HK1 | Cycle index | ATP | ADP | Non-conv. |
|---|---|---|---|---|---|---|---|---|
| r 0, 5 µM | 7376 | 0 | 0% | 3717 | 0.063 | 7472 | 28.8 | 0 |
| r 0.5, 5 µM | 4038 | 3333 | 45.2% | 2048 | 0.074 | 7445 | 56.3 | 0 |
| r 0.5, 1 µM | 3799 | 3569 | 48.4% | 1928 | 0.095 | 7395 | 105.2 | 1 |
| r 0, 1 µM | 7370 | 0 | 0% | 3714 | 0.080 | 7428 | 73.2 | 1 |
| *carbon-free r 0.5, 5 µM (`24cb571`)* | 3942 | 3431 | 46.5% | 2000 | 0.066 | 7462 | 39.1 | 0 |
| *carbon-free r 0.5, 1 µM (`24cb571`)* | 3621 | 3750 | 50.9% | 1839 | 0.076 | 7437 | 63.8 | 0 |

- **Kober's competition is reproduced.** Respiration takes about half of the demand-set ATP
  flux, and glycolysis and glucose uptake fall roughly one-for-one. Lowering F26BP from 5 to
  1 µM moves 3.2 percentage points of ATP supply from glycolysis to respiration, the direction
  Kober reports.
- **The physiological operating point moves a little.** The cycle index rises by 0.011 (5 µM)
  and 0.015 (1 µM), and ADP roughly doubles at 5 µM. Both shifts are larger than with
  carbon-free respiration.
- **Respiration does not rescue the high-NADPH collapse once it needs carbon.**
  - Without respiration, 16–17 cells per ATP level collapse at every ATP demand: ATP → 0, the
    R5P sink runs backward, and the cycle index is undefined (gray in the viewer).
  - Carbon-free respiration held these cells at ATP ≈ 7.4 mM.
  - With carbon-linked respiration they fall into a different trapped state. Taking the median
    over the 187 cells at NADPH φ ≥ 0.29 with ATP < 100 µM: ATP 0.5 µM, AMP 7.4 mM, S7P 29 mM,
    Pi 56 µM. Phosphate is locked in S7P and the adenylate pool has drained to AMP.
  - G6PD and the R5P sink stay non-negative there, so the cycle index is defined and **the
    viewer colors these cells instead of graying them out**. Check `ATP_uM` before reading a
    high-NADPH cell.
  - The carbon-free rescue was therefore an artifact of free ATP.
- **The low-demand F16BP trap expands.** At the lowest ATP demands, cells pile up about 11 mM
  F16BP with Pi in the low µM range. These are glycolysis's Harden–Young-type phosphate trap,
  which already appears on main at φ_ATP = 0.01. ATP < 1 mM cells per ATP level, from lowest
  demand up (the steady 16–17 are the high-NADPH cells):

  | ATP demand | 0.010 | 0.014 | 0.019 | 0.026 | 0.036 | ≥ 0.050 |
  |---|---|---|---|---|---|---|
  | r 0, 5 µM | 99 | 28 | 17 | 17 | 17 | 16–17 |
  | r 0.5, 5 µM | 140 | 137 | 123 | 96 | 17 | 16–18 |
