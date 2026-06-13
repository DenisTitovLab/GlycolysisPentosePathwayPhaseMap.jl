# metab_panel.jl — literature box-plot panel + dynamic model diamonds

# Free-aldehyde corrections, mirrored from PentosePhosphatePathway.jl/src/model_parameters.jl
# (GAP_correction = 1 + 1/0.0297 ≈ 34.7; DHAP_correction = 1 + 1/2.24 ≈ 1.45). Only the free
# aldehyde is kinetically active, so literature totals are divided by these to match the model's
# pool basis (as helper_files/import_data.jl does). Mirrored, not imported, to keep the explorer
# PPP-free; these are fixed thermodynamic constants with negligible drift risk.
const _GAP_CORRECTION  = 1 + 1 / 0.0297
const _DHAP_CORRECTION = 1 + 1 / 2.24

"""
    load_metab_lit(csv_path; order = _PANEL_METAB_ORDER) -> Dict{Symbol, Vector{Float64}}

Read the literature metabolite-concentration table and return, per panel metabolite, the vector
of measured values **in µM** (missing entries dropped). GAP and DHAP are divided by the
free-aldehyde corrections so they sit on the model's pool basis. A metabolite whose column is
absent or all-missing yields an empty vector. Stays in µM (no µM→M conversion) to match the
cached `<sym>_uM` pool columns the panel's model diamonds come from.
"""
function load_metab_lit(csv_path; order = _PANEL_METAB_ORDER)
    df = CSV.read(csv_path, DataFrame)
    cols = names(df)
    lit = Dict{Symbol, Vector{Float64}}()
    for sym in order
        col = string(sym)
        vals = col in cols ? Float64.(collect(skipmissing(df[!, col]))) : Float64[]
        sym === :GAP  && (vals ./= _GAP_CORRECTION)
        sym === :DHAP && (vals ./= _DHAP_CORRECTION)
        lit[sym] = vals
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
    ylims!(ax, 1e-3, 1e5)
    ax.yscale = log10
    ax.ylabel = "[metabolite], µM"
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
