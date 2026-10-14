####################### render the explorer screenshot for the README #######################
# Builds the interactive explorer figure offscreen (CairoMakie, no GL window) and saves a PNG.
# Shows the explorer's opening view: the physiological-demand cell at the default ATP level.

using CairoMakie, DataFrames, GlycolysisPentosePathwayPhaseMap
const G = GlycolysisPentosePathwayPhaseMap
CairoMakie.activate!()

df = G.load_grid(G.default_atp_grid())
fig, selected, atp_level = G.build_explorer(df)
in_, ir = selected[]
row = only(eachrow(df[(df.atpase_frac .== atp_level[]) .& (df.i_nadph .== in_) .& (df.i_r5p .== ir), :]))
println("opening cell (nadph=", in_, ", r5p=", ir, ") at ATP ", atp_level[], ": ", row.mode,
        ", cycle_index=", row.cycle_index)

out = joinpath(@__DIR__, "..", "docs", "src", "assets", "explorer.png")
mkpath(dirname(out))
save(out, fig; px_per_unit = 2)
println("wrote ", abspath(out))
