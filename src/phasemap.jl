# phasemap.jl — static Layout-A phase map renderer + public phasemap() entry

"""
    render_phase_map!(fig, df; reps) -> fig

Layout A: a cycle-index heatmap over the NADPH × R5P demand grid (left, full height) plus a
glycolysis+PPP flux-map schematic for each `reps` operating point stacked on the right. `df`
is the `run_demand_grid` output; `reps` is a vector of (label, i_nadph, i_r5p) tuples (two by
default: physiological demand + strongest reverse). Non-converged cells show as NaN (the colormap's nan color); all panels
share one flux-width scale.
"""
function render_phase_map!(fig, df; reps)
    nx = maximum(df.i_r5p); ny = maximum(df.i_nadph)
    Z = fill(NaN, nx, ny)                       # [r5p, nadph]; NaN = masked
    for r in eachrow(df)
        r.retcode == :Terminated || continue    # mask non-converged cells (spec)
        Z[r.i_r5p, r.i_nadph] = r.cycle_index
    end
    xs = sort(unique(df.r5p_phi)); ys = sort(unique(df.nadph_phi))

    # span all schematic rows so the heatmap stays full-height alongside the stacked panels
    ax = Axis(fig[1:length(reps), 1]; xscale = log10, yscale = log10,
              xlabel = "R5P demand (fraction of R5P supply, 1.5·RPI_cap)",
              ylabel = "NADPH demand (fraction of NADPH supply, 6·RPI_cap)",
              title = "Pentose-cycle index")
    hm = heatmap!(ax, log_cell_edges(xs), log_cell_edges(ys), Z; colormap = :RdBu,
                  colorrange = (-1, 1), nan_color = :gray80)
    Colorbar(fig[1:length(reps), 2], hm;
             label = "recycled fraction (1=cycle, 0=linear, <0=reverse)")
    # sign-flip contour (mode boundary) where cycle_index crosses 0
    contour!(ax, xs, ys, Z; levels = [0.0], color = :black, linewidth = 2)

    # collect each rep's fluxes first so all panels share ONE width scale: maxflux = the
    # largest |flux| over every edge of every rep.
    reps_fluxes = map(reps) do (label, i_n, i_r)
        row = first(filter(r -> r.i_nadph == i_n && r.i_r5p == i_r, eachrow(df)))
        (label, i_n, i_r, row_to_fluxes(row))
    end
    # skip NaN edges so a non-converged rep cell (masked fluxes) can't poison the shared
    # scale; a NaN maxflux would otherwise make every panel's arrow widths NaN.
    maxflux = maximum(reps_fluxes; init = 1.0) do (_, _, _, fluxes)
        maximum((abs(getfield(fluxes, k)) for k in _PPP_FLUX_KEYS); init = 0.0) do v
            isnan(v) ? 0.0 : v
        end
    end
    maxflux = (maxflux == 0 || isnan(maxflux)) ? 1.0 : maxflux

    for (k, (label, i_n, i_r, fluxes)) in enumerate(reps_fluxes)
        scatter!(ax, [xs[i_r]], [ys[i_n]]; marker = :rect, markersize = 14,
                 color = :transparent, strokecolor = :black, strokewidth = 2)
        text!(ax, xs[i_r], ys[i_n]; text = string(k), align = (:center, :center),
              fontsize = 18, font = :bold)
        axk = Axis(fig[k, 3])
        draw_ppp_flux_map!(axk, fluxes; title = "$(k): $(label)", maxflux = maxflux)
    end
    return fig
end

# Default phasemap panels: the physiological-demand cell (labelled with its mode) and the
# converged cell with the most negative cycle index (the strongest reverse operation).
function _default_reps(df)
    in_, ir = physiological_cell(df)
    phys = only(eachrow(df[(df.i_nadph .== in_) .& (df.i_r5p .== ir), :]))
    conv = filter(r -> r.retcode == :Terminated && !isnan(r.cycle_index), eachrow(df))
    rev = conv[argmin([r.cycle_index for r in conv])]
    return [("physiological demand ($(phys.mode))", in_, ir),
            ("reverse (strongest)", rev.i_nadph, rev.i_r5p)]
end

"""
    phasemap(; grid = default_demand_grid(), save = nothing,
               reps = nothing, size = (1500, 850)) -> Figure

Render the static Layout-A phase map (cycle-index heatmap + two representative flux schematics)
from a 2-D demand grid CSV (defaults to the one shipped in `data/`). Renders with CairoMakie;
writes a PNG when `save` is a path. `reps` defaults to the cell nearest physiological demand
(PPP.jl's default NADPH and R5P sinks) and the most strongly reverse converged cell.
"""
function phasemap(; grid = default_demand_grid(), save = nothing, reps = nothing,
                    size = (1500, 850))
    CairoMakie.activate!()
    df = load_grid(grid)
    reps === nothing && (reps = _default_reps(df))
    fig = Figure(; size = size)
    render_phase_map!(fig, df; reps = reps)
    save === nothing || Makie.save(save, fig)
    return fig
end
