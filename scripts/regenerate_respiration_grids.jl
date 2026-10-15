##########################################################################################
#        regenerate_respiration_grids.jl — respiration x F26BP what-if grid set          #
##########################################################################################
#
# Solves the phase-map grids on the model with an ATPSynthase respiration pseudo-enzyme
# (ext/respiration.jl) for the 2x2 respiration x F26BP comparison. The r = 0, 5 µM corner is the
# shipped data/ grid pair; this script writes the other three corners to data/respiration/<cond>/,
# each with its own provenance.toml ([grids] only).
#
# Step 1 is a drift check: the r = 0, 5 µM grids solved through the respiration code path must
# reproduce the shipped data/ grids on every shared column, or nothing is written.
#
# Run from the PPP_Experiments environment (it `develop`s this package and PPP.jl by path):
#   julia --project=/home/james/projects/PPP_Experiments scripts/regenerate_respiration_grids.jl

using GlycolysisPentosePathwayPhaseMap, PentosePhosphatePathway
using CSV, DataFrames, Printf
const G = GlycolysisPentosePathwayPhaseMap

const RESP_DIR   = joinpath(G.DATA_DIR, "respiration")
const CONDITIONS = [(name = "r0.5_f26bp5uM", resp_frac = 0.5, f26bp_uM = 5.0),
                    (name = "r0.5_f26bp1uM", resp_frac = 0.5, f26bp_uM = 1.0),
                    (name = "r0_f26bp1uM",   resp_frac = 0.0, f26bp_uM = 1.0)]

##########################################################################################
#                                   1. Drift check                                       #
##########################################################################################

# Largest relative difference over the columns two grids share; Inf if a non-numeric column
# or a NaN pattern differs.
function max_rel_diff(new, old)
    worst = 0.0
    for c in intersect(names(new), names(old))
        a, b = new[!, c], old[!, c]
        if eltype(a) <: Real && eltype(b) <: Real
            isequal(isnan.(a), isnan.(b)) || return Inf
            for (x, y) in zip(a, b)
                isnan(x) && continue
                worst = max(worst, x == y ? 0.0 : abs(x - y) / max(abs(x), abs(y)))
            end
        else
            isequal(string.(a), string.(b)) || return Inf
        end
    end
    return worst
end

drift_dir = mktempdir()
drift = regenerate_grid(; atp_out = joinpath(drift_dir, "atp_demand_grid.csv"),
                          demand_out = joinpath(drift_dir, "demand_grid.csv"),
                          resp_frac = 0.0, f26bp_uM = 5.0)
for (new, old) in ((drift.atp_out, G.default_atp_grid()), (drift.demand_out, G.default_demand_grid()))
    d = max_rel_diff(CSV.read(new, DataFrame), CSV.read(old, DataFrame))
    @printf("drift check %-22s max rel diff = %.3g\n", basename(old), d)
    d <= 1e-9 || error("r = 0 grids drift from the shipped $(basename(old)) (max rel diff $d); " *
                       "not writing the respiration grids")
end

##########################################################################################
#                                 2. Condition grids                                     #
##########################################################################################

for c in CONDITIONS
    dir = joinpath(RESP_DIR, c.name)
    regenerate_grid(; atp_out = joinpath(dir, "atp_demand_grid.csv"),
                      demand_out = joinpath(dir, "demand_grid.csv"),
                      resp_frac = c.resp_frac, f26bp_uM = c.f26bp_uM)
end

##########################################################################################
#                                     3. Summary                                         #
##########################################################################################

# Physiological cell of each corner's 2-D grid (ATP demand 0.25 of 2·HK1_cap).
function summarize(label, atp_path, demand_path)
    atp = G.load_grid(atp_path); dem = G.load_grid(demand_path)
    in_, ir = G.physiological_cell(dem)
    row = only(dem[(dem.i_nadph .== in_) .& (dem.i_r5p .== ir), :])
    resp = hasproperty(row, :V_ATPSynthase_uM_min) ? row.V_ATPSynthase_uM_min : 0.0
    glyc = hasproperty(row, :V_PKM2_uM_min) ? row.V_PKM2_uM_min : NaN
    @printf("%-15s | nonconv 3-D %3d, 2-D %3d | undetermined 3-D %4d | glyc %7.1f resp %7.1f µM/min (resp %4.1f%%) | HK1 %6.1f | cycle %6.3f | ATP %6.0f ADP %6.1f µM\n",
            label, count(!=(:Terminated), atp.retcode), count(!=(:Terminated), dem.retcode),
            count(==("undetermined"), string.(atp.mode)), glyc, resp, 100 * resp / (glyc + resp),
            row.V_HK1_uM_min, row.cycle_index, row.ATP_uM, row.ADP_uM)
end

println()
summarize("r0_f26bp5uM", drift.atp_out, drift.demand_out)
for c in CONDITIONS
    dir = joinpath(RESP_DIR, c.name)
    summarize(c.name, joinpath(dir, "atp_demand_grid.csv"), joinpath(dir, "demand_grid.csv"))
end
