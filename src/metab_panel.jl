# metab_panel.jl — metabolite panel: free ≈ total axis + free-pool axis, dynamic model marks

# The literature table is exported by PPP_Experiments (scripts/phasemap_export/) already on the
# model's pool basis: µM of cytosolic water, GAP/DHAP as the free aldehyde, and free-pool data
# only for the _FREE_POOL_PANEL cofactors. This package does no basis conversion. The split
# mirrors PPP_Experiments' basic panel (scripts/helper_files/plot_metab_split.jl).

"Metabolite-panel y-axis label; states the basis both the data and the model marks are on."
const _METAB_YLABEL = "[metabolite], µM (cytosolic water)"

"Metabolite-panel y-limits, µM. 1e-4 keeps every non-collapsed grid cell's pools on scale."
const _METAB_YLIMS = (1e-4, 1e5)

"Allowed `method` values in literature_pools.csv."
const _METAB_METHODS = (:metabolomics, :sensor, :equilibrium)

"Extra x-room right of the last column, for an off-scale value label."
const _XPAD_RIGHT = 0.3

"One literature value: µM, whether it was derived (matched pair), and how it was measured."
const _LitPoint = @NamedTuple{value::Float64, derived::Bool, method::Symbol}

"""
    load_metab_lit(csv_path; order = _PANEL_METAB_ORDER) -> Dict{Symbol, Vector{_LitPoint}}

Read the long-format literature table (`metabolite, value_uM, reference, cell_type, derived,
method`) and return, per panel metabolite, its values as `_LitPoint`s. A panel metabolite with
no rows maps to an empty vector; rows outside `order` are ignored. Throws `ArgumentError` on a
`method` outside `_METAB_METHODS`, on a non-metabolomics row for a `_FREE_EQ_TOTAL_PANEL`
metabolite, and on a metabolomics (total-pool) row for a `_FREE_POOL_PANEL` metabolite.
"""
function load_metab_lit(csv_path; order = _PANEL_METAB_ORDER)
    df = CSV.read(csv_path, DataFrame; types = Dict(:metabolite => String, :value_uM => Float64,
                                                    :derived => Bool, :method => String))
    lit = Dict{Symbol, Vector{_LitPoint}}(sym => _LitPoint[] for sym in order)
    for r in eachrow(df)
        sym, m = Symbol(r.metabolite), Symbol(r.method)
        m in _METAB_METHODS || throw(ArgumentError(
            "literature row for $sym has method :$m; expected one of $(_METAB_METHODS)"))
        haskey(lit, sym) || continue
        if sym in _FREE_POOL_PANEL
            m === :metabolomics && throw(ArgumentError(
                "free-pool metabolite $sym has a metabolomics (total-pool) row"))
        elseif m !== :metabolomics
            throw(ArgumentError("free ≈ total metabolite $sym has a non-metabolomics row (:$m)"))
        end
        push!(lit[sym], (value = r.value_uM, derived = r.derived, method = m))
    end
    return lit
end

"""
    _free_pool_point_style(derived, method) -> NamedTuple

Marker for one free-pool data point: filled = measured, open = inferred (by us or the original
authors); circle = sensor-based, diamond = equilibrium-based. Copied from PPP_Experiments
scripts/helper_files/plot_metab_split.jl (`free_pool_point_style`) so both panels share one key.
"""
function _free_pool_point_style(derived::Bool, method::Symbol)
    method === :equilibrium && return (label = "equilibrium-based", marker = :diamond,
        color = :transparent, strokecolor = :black, strokewidth = 1.5)
    derived && return (label = "derived free (sensor ratio)", marker = :circle,
        color = :transparent, strokecolor = :black, strokewidth = 1.5)
    return (label = "measured free (sensor)", marker = :circle, color = :black,
        strokecolor = :black, strokewidth = 0.0)
end

"""
    _offscale_marks(pools; lims = _METAB_YLIMS) -> (; on, below, above, labels, label_pos)

Split one axis's model pools (µM, column i = position i) into in-range diamond positions and
floor/ceiling arrow positions with value labels. Zero and negative solver values count as below
the floor; non-finite values (solver-error rows) draw nothing.
"""
function _offscale_marks(pools; lims = _METAB_YLIMS)
    lo, hi = lims
    on, below, above, label_pos = Point2f[], Point2f[], Point2f[], Point2f[]
    labels = String[]
    for (i, v) in enumerate(pools)
        isfinite(v) || continue
        if lo <= v <= hi
            push!(on, Point2f(i, v)); continue
        end
        y = v < lo ? lo * 1.5 : hi / 1.5
        push!(v < lo ? below : above, Point2f(i, y))
        push!(labels, string(round(v; sigdigits = 2))); push!(label_pos, Point2f(i + 0.12, y))
    end
    return (; on, below, above, labels, label_pos)
end

# One metabolite axis at gl[1, col]: categorical x pinned to the columns (sparse reactive marks
# must never auto-limit it), log y with the limits set BEFORE the log switch.
function _metab_axis(gl, col, title, labels)
    ax = Axis(gl[1, col]; title, xticks = (collect(eachindex(labels)), labels),
              xticklabelrotation = pi / 3)
    ylims!(ax, _METAB_YLIMS...)
    ax.yscale = log10
    xlims!(ax, 0.5, length(labels) + 0.5 + _XPAD_RIGHT)
    return ax
end

# Model marks lifted off `pools_obs`: a fixed set of four plots whose data change on reselection.
function _draw_model_marks!(ax, pools_obs)
    marks = lift(_offscale_marks, pools_obs)
    scatter!(ax, lift(m -> m.on, marks); marker = :diamond, markersize = 21, color = :purple,
             strokecolor = :black, strokewidth = 0.5)
    scatter!(ax, lift(m -> m.below, marks); marker = :dtriangle, markersize = 14, color = :purple)
    scatter!(ax, lift(m -> m.above, marks); marker = :utriangle, markersize = 14, color = :purple)
    text!(ax, lift(m -> m.label_pos, marks); text = lift(m -> m.labels, marks), fontsize = 10,
          color = :purple, align = (:left, :center))
end

# Static metabolomics boxes, one boxplot! call grouped by integer x-position.
function _draw_lit_boxes!(ax, lit, metabs)
    box_x = Float64[]; box_y = Float64[]
    for (i, sym) in enumerate(metabs), p in lit[sym]
        push!(box_x, i); push!(box_y, p.value)
    end
    isempty(box_x) || boxplot!(ax, box_x, box_y; width = 0.6, color = (:gray70, 0.6),
                               strokecolor = :gray30, strokewidth = 1, whiskerwidth = 0.6,
                               markersize = 4)
end

# Static free-pool points, spread over ±0.2 of the column when there are several.
function _draw_free_pool_points!(ax, lit, metabs)
    for (x, sym) in enumerate(metabs)
        pts = lit[sym]
        isempty(pts) && continue
        xs = length(pts) == 1 ? [Float64(x)] : collect(range(x - 0.2, x + 0.2, length = length(pts)))
        for (xi, p) in zip(xs, pts)
            s = _free_pool_point_style(p.derived, p.method)
            scatter!(ax, [xi], [p.value]; marker = s.marker, color = s.color,
                     strokecolor = s.strokecolor, strokewidth = s.strokewidth, markersize = 11)
        end
    end
end

function _metab_legend!(pos)
    styles = [_free_pool_point_style(false, :sensor), _free_pool_point_style(true, :sensor),
              _free_pool_point_style(false, :equilibrium)]
    elems = [PolyElement(color = (:gray70, 0.6), strokecolor = :gray30),
             MarkerElement(marker = :diamond, color = :purple, strokecolor = :black,
                           strokewidth = 0.5, markersize = 16),
             MarkerElement(marker = :dtriangle, color = :purple, markersize = 12),
             (MarkerElement(marker = s.marker, color = s.color, strokecolor = s.strokecolor,
                            strokewidth = s.strokewidth, markersize = 11) for s in styles)...]
    labels = ["literature (metabolomics)", "selected cell", "selected cell, off scale",
              (s.label for s in styles)...]
    return Legend(pos, elems, labels; framevisible = false, labelsize = 11, tellheight = false)
end

"""
    draw_metab_split!(gl, lit, metab_obs) -> (; ax_left, ax_right)

Render the metabolite panel into the `GridLayout` `gl`: left axis `_FREE_EQ_TOTAL_PANEL`
(static metabolomics boxes), right axis `_FREE_POOL_PANEL` (static free-pool points; a column
with none reads "no data"), legend in `gl[1, 3]` outside both data areas. `metab_obs` is a
length-`_PANEL_METAB_ORDER` Observable of the selected cell's pools (µM); its purple diamonds,
and labelled floor/ceiling triangles for off-scale values, update on reselection without
adding or removing plots. Linked log y over `_METAB_YLIMS`; widths follow metabolite counts.
"""
function draw_metab_split!(gl, lit, metab_obs)
    left, right = collect(_FREE_EQ_TOTAL_PANEL), collect(_FREE_POOL_PANEL)
    nl = length(left)
    right_labels = [isempty(lit[s]) ? "$s\nno data" : string(s) for s in right]
    ax_left = _metab_axis(gl, 1, "free ≈ total (extraction metabolomics)", string.(left))
    ax_left.ylabel = _METAB_YLABEL
    ax_right = _metab_axis(gl, 2, "free pool (free-pool data only)", right_labels)
    linkyaxes!(ax_left, ax_right)
    hideydecorations!(ax_right; grid = false, ticks = false)
    colsize!(gl, 1, Auto(length(left))); colsize!(gl, 2, Auto(length(right))); colgap!(gl, 12)
    _draw_lit_boxes!(ax_left, lit, left)
    _draw_free_pool_points!(ax_right, lit, right)
    _draw_model_marks!(ax_left, lift(p -> p[1:nl], metab_obs))
    _draw_model_marks!(ax_right, lift(p -> p[nl+1:end], metab_obs))
    _metab_legend!(gl[1, 3])
    return (; ax_left, ax_right)
end
