module GlycolysisPentosePathwayPhaseMap

using Makie, CSV, DataFrames
import GLMakie, CairoMakie   # backends activated per-entry; NOT `using` (avoids dual re-export clash)

export explore, phasemap, regenerate_grid

const DATA_DIR = joinpath(@__DIR__, "..", "data")
default_atp_grid()    = joinpath(DATA_DIR, "atp_demand_grid.csv")
default_demand_grid() = joinpath(DATA_DIR, "demand_grid.csv")
metab_lit_path()      = joinpath(DATA_DIR, "Metabolite_Concentrations.csv")

include("grid_io.jl")
include("flux_map.jl")
include("metab_panel.jl")
include("phasemap.jl")
include("explorer.jl")

"""
    regenerate_grid(; kwargs...)

Recompute the shipped grid CSVs by solving the kinetic ODE model. Only available when the
recompute extension is active — load it first with `using PentosePhosphatePathway`.
"""
function regenerate_grid(; kwargs...)
    error("regenerate_grid needs the recompute extension. Load the kinetic model first:\n" *
          "    using PentosePhosphatePathway\n" *
          "(this pulls in the solver stack and activates PentosePhosphatePathwayExt).")
end

end # module
