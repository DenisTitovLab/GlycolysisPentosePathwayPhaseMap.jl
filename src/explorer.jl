# explorer.jl — interactive demand-grid flux explorer (GLMakie window) + offscreen builder

# The heatmap draws cell (i_r5p, i_nadph) centred at that index pair, so a click snaps to the
# nearest integer position, clamped onto the grid. Returns (i_r5p, i_nadph).
click_to_cell(nx, ny, pos) =
    (clamp(round(Int, pos[1]), 1, nx), clamp(round(Int, pos[2]), 1, ny))

# O(1) row lookup keyed (atpase_frac, i_nadph, i_r5p), built once. Exact == on the CSV-parsed
# Float64 atpase_frac is safe: the slider values come from the same parsed column.
_build_row_index(df) = Dict((r.atpase_frac, r.i_nadph, r.i_r5p) => r for r in eachrow(df))

# Per-ATP-level cycle-index matrices [i_r5p, i_nadph], built once. Masked/non-converged cells
# stay NaN (renderers mask on the :Terminated retcode). Replaces a full df rescan per slider tick.
function _build_zdict(df, nx, ny)
    zd = Dict(atp => fill(NaN, nx, ny) for atp in unique(df.atpase_frac))
    for r in eachrow(df)
        r.retcode == :Terminated || continue
        zd[r.atpase_frac][r.i_r5p, r.i_nadph] = r.cycle_index
    end
    return zd
end

"""
    build_explorer(df) -> (fig, selected, atp_level)

Build the explorer over the 3-D grid `df` (`retcode` already Symbol-normalized). The grid must carry
an `atpase_frac` column and a `<sym>_uM` column for every metabolite in `_PANEL_METAB_ORDER`;
otherwise it throws an `ArgumentError` saying to regenerate the grid. An ATP-demand slider re-slices the cached grid: the heatmap Z + cycle-index
contour lift off the selected ATP level, and the network redraws on EITHER an ATP-level change or a
cell click. Returns the Figure, the `selected` Observable (keyed `(i_nadph, i_r5p)`), and the
`atp_level` Observable. Does NOT call `display` (offscreen-testable).
"""
function build_explorer(df)
    _require_columns(df, ["atpase_frac"],
                     "the explorer needs the 3-D ATP grid (atp_demand_grid.csv)")
    _require_columns(df, ["$(s)_uM" for s in _PANEL_METAB_ORDER],
                     "the grid predates the cofactor pool columns; regenerate it with regenerate_grid()")
    # ATP levels and demand axes come straight from the CSV-parsed df, so slider values and the
    # filter compare bit-identical Float64s (exact `==` is safe).
    atp_levels = sort(unique(df.atpase_frac))
    xs = sort(unique(df.r5p_phi)); ys = sort(unique(df.nadph_phi))
    nx = maximum(df.i_r5p); ny = maximum(df.i_nadph)

    row_index = _build_row_index(df)
    zdict = _build_zdict(df, nx, ny)

    fig = Figure(size = (1500, 1200))
    ax_heat = Axis(fig[1, 1]; _demand_axis_kwargs(xs, ys)...,
                   title = "Pentose-cycle index — click a cell")

    # ATP-demand slider (log-spaced levels straight from the data); default = the level nearest
    # the physiological operating point's ATPase demand.
    default_atp = atp_levels[argmin(abs.(log10.(atp_levels) .-
                                         log10(_PHYSIOLOGICAL_DEMAND.atpase_frac)))]
    sg = SliderGrid(fig[3, 1:3],
                    (label = "ATP demand (% of max supply)", range = atp_levels,
                     startvalue = default_atp,
                     format = x -> string(round(x * 100, sigdigits = 2), "%")))
    atp_level = sg.sliders[1].value

    Z_obs = lift(atp -> zdict[atp], atp_level)
    # one equal-sized cell per demand level, centred on its index (see _demand_axis_kwargs)
    hm = heatmap!(ax_heat, 1:nx, 1:ny, Z_obs; colormap = :RdBu,
                  colorrange = (-1, 1), nan_color = :gray80)
    Colorbar(fig[1, 2], hm; label = "recycled fraction (1=cycle, 0=linear, <0=reverse)")
    contour!(ax_heat, 1:nx, 1:ny, Z_obs; levels = [0.0], color = :black, linewidth = 2)
    for interaction in (:rectanglezoom, :dragpan, :scrollzoom)
        deregister_interaction!(ax_heat, interaction)
    end

    ax_net = Axis(fig[1, 3])
    footer = Observable("")
    Label(fig[2, 1:3], footer; halign = :left, fontsize = 19, tellwidth = false)
    colsize!(fig.layout, 3, Relative(0.60))
    rowsize!(fig.layout, 2, Relative(0.04))   # thin footer line

    # Footer row 2, ATP slider row 3, full-width metabolite panel row 4: two linked axes
    # (free ≈ total | free pool) plus their legend, in a nested GridLayout.
    metab_gl = GridLayout(fig[4, 1:3])
    rowsize!(fig.layout, 4, Relative(0.30))
    lit = load_metab_lit(metab_lit_path())
    metab_obs = Observable(fill(NaN, length(_PANEL_METAB_ORDER)))
    metab_axes = draw_metab_split!(metab_gl, lit, metab_obs)
    for ax in metab_axes, interaction in (:rectanglezoom, :dragpan, :scrollzoom)
        deregister_interaction!(ax, interaction)
    end

    selected = Observable(physiological_cell(df))   # open on the physiological-demand cell
    marker_pos = lift(selected) do (in_, ir)
        Point2f(ir, in_)
    end
    scatter!(ax_heat, marker_pos; marker = :rect, markersize = 18,
             color = :transparent, strokecolor = :black, strokewidth = 3)
    # the exact physiological demand, which falls between grid cells; stays put while you explore
    scatter!(ax_heat, [index_position(xs, _PHYSIOLOGICAL_DEMAND.r5p_phi)],
             [index_position(ys, _PHYSIOLOGICAL_DEMAND.nadph_phi)];
             marker = :circle, markersize = 7, color = :black, strokecolor = :white,
             strokewidth = 1)

    # build the network ONCE against Observables; reselect just pushes new values (no rebuild)
    fluxes_obs = Observable(row_to_fluxes(row_index[(atp_level[], selected[][1], selected[][2])]))
    title_obs = Observable("")
    draw_ppp_flux_map!(ax_net, fluxes_obs; title = title_obs)

    # redraw-free update of network + footer + metab on EITHER an ATP-level change or a cell click
    onany(atp_level, selected) do atp, (in_, ir)
        row = row_index[(atp, in_, ir)]
        atp_pct = round(atp * 100, sigdigits = 2)
        title_obs[] =
            "ATP $(atp_pct)%  |  NADPH demand $(round(row.nadph_phi * 100, sigdigits = 2))% / " *
            "R5P demand $(round(row.r5p_phi * 100, sigdigits = 2))% of supply"
        fluxes_obs[] = row_to_fluxes(row)
        mx = maximum((abs(getfield(fluxes_obs[], k)) for k in _PPP_FLUX_KEYS); init = 0.0) do v
            isnan(v) ? 0.0 : v
        end
        mx = (mx == 0 || isnan(mx)) ? 1.0 : mx
        footer[] = row.retcode == :Terminated ?
            "ATP $(atp_pct)%  |  cell (nadph=$in_, r5p=$ir)  |  mode=$(row.mode)  |  " *
            "cycle index=$(round(row.cycle_index, sigdigits = 3))  |  " *
            "max |flux| = $(round(mx, sigdigits = 3)) µM/min" :
            "ATP $(atp_pct)%  |  cell (nadph=$in_, r5p=$ir)  |  NON-CONVERGED ($(row.retcode))"
        metab_obs[] = row_to_pools(row)
    end

    on(events(fig.scene).mousebutton) do ev
        if ev.button == Mouse.left && ev.action == Mouse.press && is_mouseinside(ax_heat.scene)
            ir, in_ = click_to_cell(nx, ny, mouseposition(ax_heat.scene))
            selected[] = (in_, ir)
        end
    end

    notify(selected)                               # draw initial cell at the default ATP level
    return fig, selected, atp_level
end

"""
    explore(; grid = default_atp_grid())
    explore(grid::AbstractString)

Open the interactive GLMakie explorer over a cached 3-D ATP×NADPH×R5P grid CSV (defaults to the
one shipped in `data/`). Click heatmap cells and drag the ATP-demand slider to re-slice the cached
grid — no ODE solve. Blocks until the window is closed.
"""
function explore(; grid = default_atp_grid())
    GLMakie.activate!()
    df = load_grid(grid)
    fig, _, _ = build_explorer(df)
    screen = GLMakie.display(fig)
    @info "Explorer open — click heatmap cells, drag the ATP slider. Close the window to exit."
    wait(screen)
    return nothing
end
explore(grid::AbstractString) = explore(; grid = grid)
