# explorer.jl — interactive demand-grid flux explorer (GLMakie window) + offscreen builder

# `mouseposition` on a log10-scaled axis returns coordinates in the axis's *transformed*
# (log10) space, not raw demand values. Convert back to data space before `nearest_cell`
# (which applies log10 itself). Returns (i_r5p, i_nadph).
logclick_to_cell(xs, ys, pos) = nearest_cell(xs, ys, 10.0^pos[1], 10.0^pos[2])

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

Build the explorer over the 3-D grid `df` (must carry an `atpase_frac` column; `retcode` already
Symbol-normalized). An ATP-demand slider re-slices the cached grid: the heatmap Z + cycle-index
contour lift off the selected ATP level, and the network redraws on EITHER an ATP-level change or a
cell click. Returns the Figure, the `selected` Observable (keyed `(i_nadph, i_r5p)`), and the
`atp_level` Observable. Does NOT call `display` (offscreen-testable).
"""
function build_explorer(df)
    # ATP levels and demand axes come straight from the CSV-parsed df, so slider values and the
    # filter compare bit-identical Float64s (exact `==` is safe).
    atp_levels = sort(unique(df.atpase_frac))
    xs = sort(unique(df.r5p_phi)); ys = sort(unique(df.nadph_phi))
    nx = maximum(df.i_r5p); ny = maximum(df.i_nadph)

    # cycle-index Z matrix [r5p, nadph] for one ATP level; masked/non-converged cells stay NaN
    zmatrix(atp) = begin
        Z = fill(NaN, nx, ny)
        for r in eachrow(df)
            (r.atpase_frac == atp && r.retcode == :Terminated) || continue
            Z[r.i_r5p, r.i_nadph] = r.cycle_index
        end
        Z
    end

    fig = Figure(size = (1500, 1200))
    ax_heat = Axis(fig[1, 1]; xscale = log10, yscale = log10,
                   xlabel = "R5P demand (fraction of R5P supply, 1.5·RPI_cap)",
                   ylabel = "NADPH demand (fraction of NADPH supply, 6·RPI_cap)",
                   title = "Pentose-cycle index — click a cell")

    # ATP-demand slider (log-spaced levels straight from the data); default = nearest 10%.
    default_atp = atp_levels[argmin(abs.(log10.(atp_levels) .- log10(0.10)))]
    sg = SliderGrid(fig[3, 1:3],
                    (label = "ATP demand (% of max supply)", range = atp_levels,
                     startvalue = default_atp,
                     format = x -> string(round(x * 100, sigdigits = 2), "%")))
    atp_level = sg.sliders[1].value

    Z_obs = lift(zmatrix, atp_level)
    # log-space cell edges (not the length-n centers) so cells render equal-sized on the log axis
    hm = heatmap!(ax_heat, log_cell_edges(xs), log_cell_edges(ys), Z_obs; colormap = :RdBu,
                  colorrange = (-1, 1), nan_color = :gray80)
    Colorbar(fig[1, 2], hm; label = "recycled fraction (1=cycle, 0=linear, <0=reverse)")
    contour!(ax_heat, xs, ys, Z_obs; levels = [0.0], color = :black, linewidth = 2)
    for interaction in (:rectanglezoom, :dragpan, :scrollzoom)
        deregister_interaction!(ax_heat, interaction)
    end

    ax_net = Axis(fig[1, 3])
    footer = Observable("")
    Label(fig[2, 1:3], footer; halign = :left, fontsize = 19, tellwidth = false)
    colsize!(fig.layout, 3, Relative(0.60))
    rowsize!(fig.layout, 2, Relative(0.04))   # thin footer line

    # Footer row 2, ATP slider row 3, NEW full-width metabolite panel row 4 (slider above it).
    ax_metab = Axis(fig[4, 1:3])
    rowsize!(fig.layout, 4, Relative(0.30))
    lit = load_metab_lit(metab_lit_path())
    metab_obs = Observable(fill(NaN, length(_PANEL_METAB_ORDER)))
    draw_metab_panel!(ax_metab, lit, metab_obs)
    for interaction in (:rectanglezoom, :dragpan, :scrollzoom)
        deregister_interaction!(ax_metab, interaction)
    end

    selected = Observable((cld(ny, 2), cld(nx, 2)))  # start near the middle of the demand grid
    marker_pos = lift(selected) do (in_, ir)
        Point2f(xs[ir], ys[in_])
    end
    scatter!(ax_heat, marker_pos; marker = :rect, markersize = 18,
             color = :transparent, strokecolor = :black, strokewidth = 3)

    # redraw network + footer on EITHER an ATP-level change or a cell click
    onany(atp_level, selected) do atp, (in_, ir)
        row = first(filter(r -> r.atpase_frac == atp && r.i_nadph == in_ && r.i_r5p == ir,
                           eachrow(df)))
        fluxes = row_to_fluxes(row)
        mx = maximum((abs(getfield(fluxes, k)) for k in _PPP_FLUX_KEYS); init = 0.0) do v
            isnan(v) ? 0.0 : v
        end
        mx = (mx == 0 || isnan(mx)) ? 1.0 : mx
        empty!(ax_net)
        atp_pct = round(atp * 100, sigdigits = 2)
        ttl = "ATP $(atp_pct)%  |  NADPH demand $(round(row.nadph_phi * 100, sigdigits = 2))% / " *
              "R5P demand $(round(row.r5p_phi * 100, sigdigits = 2))% of supply"
        # no maxflux override: draw_ppp_flux_map! self-normalizes glycolysis and PPP on separate
        # width scales (PPP max width = ½ glycolysis), so the PPP branch stays legible per cell.
        draw_ppp_flux_map!(ax_net, fluxes; title = ttl)
        footer[] = row.retcode == :Terminated ?
            "ATP $(atp_pct)%  |  cell (nadph=$in_, r5p=$ir)  |  mode=$(row.mode)  |  " *
            "cycle index=$(round(row.cycle_index, sigdigits = 3))  |  " *
            "max |flux| = $(round(mx, sigdigits = 3)) µM/min" :
            "ATP $(atp_pct)%  |  cell (nadph=$in_, r5p=$ir)  |  NON-CONVERGED ($(row.retcode))"
        metab_obs[] = row_to_pools(row)
    end

    on(events(fig.scene).mousebutton) do ev
        if ev.button == Mouse.left && ev.action == Mouse.press && is_mouseinside(ax_heat.scene)
            ir, in_ = logclick_to_cell(xs, ys, mouseposition(ax_heat.scene))
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
