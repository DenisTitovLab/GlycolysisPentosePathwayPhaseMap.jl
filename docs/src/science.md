# The science

This page explains the demand convention, the cycle-index readout, the mode classification, and
the biological context the phase map recapitulates.

## Supply/demand convention

The phase map is built on the integrated glycolysis + PPP kinetic model (the integrated model is
required because it *closes* the pentose cycle: the F6P and GAP produced by the non-oxidative PPP
re-enter glycolysis and regenerate G6P). Three demand levers drive the system, each implemented as
a tunable maximal-rate sink:

- **NADPH demand** — redox / biosynthetic draw on NADPH,
- **R5P demand** — ribose / nucleotide draw (a PRPS1-anchored sink),
- **ATP demand** — energy draw / glycolytic pull (held as a fixed background per run, and varied
  by the explorer's ATP slider).

Each lever is set as a fraction **φ** of that species' *sustainable supply ceiling* `S_X`, so all
three axes mean the same thing — "percent of what the pathway can actually make." The ceilings are
anchored on the carbon bottleneck, ribose-5-phosphate isomerase (RPI):

| Species | Supply ceiling `S_X` |
|---|---|
| NADPH | `6 · RPI_cap`  (= `2 · V_ox,max`) |
| R5P   | `1.5 · RPI_cap` |
| ATP   | `2 · HK1_cap`  |

with `RPI_cap = RPI_Vmax · RPI_Conc`. RPI is the bottleneck because it carries one third of the
ribulose-5-phosphate in the full cycle, so the maximum *sustainable* oxidative flux is
`3 · RPI_cap` and the NADPH ceiling is `2 · (3 · RPI_cap) = 6 · RPI_cap` (the oxidative enzymes
are far over-provisioned relative to RPI/RPE clearance, so they are not the operative ceiling).
R5P has two sources — RPI forward and the reverse non-oxidative ribose route — but the X5P
co-produced by the reverse route must return through RPE → Ru5P → RPI, so two thirds of R5P still
funnels through RPI and the sustainable R5P ceiling is `1.5 · RPI_cap`. ATP demand anchors on the
glycolytic ceiling `2 · HK1_cap`.

The reversible non-oxidative reactions (RPI, RPE, the two transketolase reactions, and
transaldolase) **change sign** between regimes; their net direction is the physical readout of
linear ↔ pentose-cycle ↔ reverse PPP operation.

## The pentose-cycle index

The single scalar summarizing each operating point is the **pentose-cycle index**

```
cycle index = (V_ox − V_R5Pase) / V_ox
```

where `V_ox` is the oxidative flux through G6PD and `V_R5Pase` is the R5P-export (demand) flux.
It is the fraction of oxidative-PPP-generated pentose carbon that is **recycled** back into
glycolysis rather than exported as ribose:

- **≈ 1** — full pentose cycle: oxidative carbon is recycled (the pathway runs to make NADPH),
- **≈ 0** — linear PPP: oxidative carbon leaves as R5P,
- **< 0** — reverse PPP: R5P is produced non-oxidatively, faster than oxidative supply,
- **NaN** — undefined when `V_ox ≈ 0`.

## Mode classification

Each cell is classified into one of four modes from its fluxes:

- **`:cycle`** — forward pentose cycle (index ≈ 1),
- **`:linear`** — linear oxidative-then-export operation (index ≈ 0),
- **`:reverse`** — reverse / non-oxidative R5P production (index < 0),
- **`:undetermined`** — oxidative flux too small to classify (index undefined).

Non-converged operating points are masked in the heatmaps (and counted), never silently dropped.

## Biological context

The map recapitulates the phenomenology reported by Feng et al., *"Nonoxidative pentose phosphate
pathway regulates CD8⁺ T cell immunity by maintaining NADPH homeostasis,"* **PNAS** 2026,
123(8):e2526325123: as demand shifts, the model spontaneously switches the PPP between linear,
forward-cycle, and reverse operation, and the non-oxidative branch carries the load of maintaining
NADPH homeostasis and ribose supply.
