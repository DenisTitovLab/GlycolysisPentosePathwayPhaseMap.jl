# metab_panel.jl — literature box-plot panel + dynamic model diamonds

# The literature table is exported by PPP_Experiments (scripts/phasemap_export/) already on the
# model's pool basis: µM of cytosolic water, free NADP(H) (NADP derived as free NADPH / R),
# GAP/DHAP as the free aldehyde. This package does no basis conversion.

"Metabolite-panel y-axis label; states the basis both the boxes and the diamonds are on."
const _METAB_YLABEL = "[metabolite], µM (cytosolic water; NADP(H) free)"

"Metabolite-panel y-limits, µM. 1e-4 keeps every non-collapsed grid cell's pools on scale."
const _METAB_YLIMS = (1e-4, 1e5)

"""
    load_metab_lit(csv_path; order = _PANEL_METAB_ORDER) -> Dict{Symbol, Vector{Float64}}

Read the long-format literature table (`metabolite, value_uM, reference, cell_type, derived`)
and return, per panel metabolite, its values in µM. Values are used as shipped: the table is
already on the model's pool basis. A panel metabolite with no rows maps to an empty vector;
rows for metabolites outside `order` are ignored.
"""
function load_metab_lit(csv_path; order = _PANEL_METAB_ORDER)
    df = CSV.read(csv_path, DataFrame; types = Dict(:metabolite => String, :value_uM => Float64))
    lit = Dict{Symbol, Vector{Float64}}(sym => Float64[] for sym in order)
    for r in eachrow(df)
        sym = Symbol(r.metabolite)
        haskey(lit, sym) && push!(lit[sym], r.value_uM)
    end
    return lit
end

"""
    draw_metab_panel!(ax, lit, metab_obs; order = _PANEL_METAB_ORDER) -> ax

Render the metabolite panel into `ax`: a STATIC literature box-and-whisker per metabolite (drawn
once from `lit`, the Dict from `load_metab_lit`) plus a DYNAMIC purple diamond per metabolite whose
position `lift`s off `metab_obs` (a length-`order` Vector of the selected cell's pools, µM). Log y
in µM, categorical x. Backend-neutral (CairoMakie static / GLMakie interactive). Non-finite or
≤0 pools draw no diamond. Reselection only updates the diamonds; the boxes stay put.
"""
function draw_metab_panel!(ax, lit, metab_obs; order = _PANEL_METAB_ORDER)
    labels = collect(string.(order))
    # set finite positive y-limits BEFORE switching to log10, so the yscale-change
    # validation never sees the default auto-limits (which include 0.0, invalid for log10)
    ylims!(ax, _METAB_YLIMS...)
    ax.yscale = log10
    ax.ylabel = _METAB_YLABEL
    ax.title  = "Steady-state metabolites — selected cell vs literature"
    ax.xticks = (1:length(order), labels)
    ax.xticklabelrotation = pi / 3

    # static literature boxes: one boxplot! call grouped by integer x-position
    box_x = Float64[]; box_y = Float64[]
    for (i, sym) in enumerate(order), v in get(lit, sym, Float64[])
        if isfinite(v) && v > 0
            push!(box_x, i); push!(box_y, v)
        end
    end
    if !isempty(box_x)
        boxplot!(ax, box_x, box_y; width = 0.6, color = (:gray70, 0.6), strokecolor = :gray30,
                 strokewidth = 1, whiskerwidth = 0.6, markersize = 4)
    end

    # dynamic model diamonds (skip non-finite / non-positive so log-y never errors)
    pts = lift(metab_obs) do pools
        Point2f[Point2f(i, pools[i]) for i in eachindex(pools) if isfinite(pools[i]) && pools[i] > 0]
    end
    scatter!(ax, pts; marker = :diamond, markersize = 21, color = :purple,
             strokecolor = :black, strokewidth = 0.5)

    axislegend(ax,
        [PolyElement(color = (:gray70, 0.6), strokecolor = :gray30),
         MarkerElement(marker = :diamond, color = :purple, markersize = 18)],
        ["literature", "selected cell"];
        position = :rt, labelsize = 11, framevisible = false)
    return ax
end
