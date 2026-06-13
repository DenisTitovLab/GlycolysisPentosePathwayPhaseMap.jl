####################### render the explorer screenshot for the README #######################
# Builds the interactive explorer figure offscreen (CairoMakie, no GL window) and saves a PNG.
# Picks a converged cell so the flux network panel shows meaningful flux.

using CairoMakie, DataFrames, GlycolysisPentosePathwayPhaseMap
const G = GlycolysisPentosePathwayPhaseMap
CairoMakie.activate!()

df = G.load_grid(G.default_atp_grid())
fig, selected, atp_level = G.build_explorer(df)

# Pick a converged cell with a strong, well-defined forward pentose cycle (cycle_index near 1,
# excluding NaN/degenerate cells) at the default ~10% ATP level, so the network redraw shows
# clear classified flux rather than a blank or undetermined cell.
atp = atp_level[]
conv = filter(r -> r.retcode == :Terminated && r.atpase_frac == atp &&
                   !isnan(r.cycle_index), eachrow(df))
best = conv[argmin([abs(r.cycle_index - 0.9) for r in conv])]
selected[] = (best.i_nadph, best.i_r5p)
println("picked cell (nadph=", best.i_nadph, ", r5p=", best.i_r5p,
        ") cycle_index=", best.cycle_index)

out = joinpath(@__DIR__, "..", "docs", "src", "assets", "explorer.png")
mkpath(dirname(out))
save(out, fig; px_per_unit = 2)
println("wrote ", abspath(out), "  (cell cycle_index=", best.cycle_index, ")")
