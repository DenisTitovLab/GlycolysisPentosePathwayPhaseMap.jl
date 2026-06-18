##########################################################################################
#                  grid_io.jl — pure flux diagnostics + grid-CSV helpers                 #
##########################################################################################

# Display order for the metabolite panel: glycolysis spine → oxidative → non-oxidative PPP
# plus the redox cofactors. Single source of truth used by row_to_pools and metab_panel.jl.
const _PANEL_METAB_ORDER = (:Glucose, :G6P, :F6P, :F16BP, :DHAP, :GAP, :ATP, :NADP, :NADPH,
                            :PGLn, :PGA, :Ru5P, :R5P, :X5P, :S7P, :E4P)

"""
    pentose_cycle_index(fluxes; tol = 1e-12) -> Float64

Fraction of oxidative-PPP-generated pentose carbon recycled rather than exported as
R5P:  (V_ox - V_R5Pase) / V_ox  with  V_ox = fluxes.V_G6PD.
≈1 full pentose cycle (max NADPH/glucose); ≈0 linear; <0 net non-oxidative R5P
production (non-ox runs reverse). Returns NaN when V_ox ≤ tol (ill-defined).
"""
function pentose_cycle_index(fluxes; tol = 1e-12)
    V_ox = fluxes.V_G6PD
    V_ox <= tol && return NaN
    return (V_ox - fluxes.V_R5Pase) / V_ox
end

"""
    classify_mode(fluxes; cycle_thresh = 0.5) -> Symbol

:cycle (index ≥ cycle_thresh), :reverse (index ≤ 0), :linear (in between), or
:undetermined (index is NaN, i.e. oxPPP flux ≈ 0).
"""
function classify_mode(fluxes; cycle_thresh = 0.5)
    ci = pentose_cycle_index(fluxes)
    isnan(ci) && return :undetermined
    ci >= cycle_thresh && return :cycle
    ci <= 0.0         && return :reverse
    return :linear
end

"""
    row_to_fluxes(row) -> NamedTuple

Build the reaction-flux NamedTuple (µM/min, signed) that `draw_ppp_flux_map!` consumes,
from one grid-CSV `DataFrameRow`. The 15 keys mirror the flux-map reactions.
"""
row_to_fluxes(row) =
    (; V_HK1 = row.V_HK1_uM_min, V_GPI = row.V_GPI_uM_min,
       V_PFKP = row.V_PFKP_uM_min, V_ALDO = row.V_ALDO_uM_min,
       V_TPI = row.V_TPI_uM_min, V_GAPDH = row.V_GAPDH_uM_min,
       V_G6PD = row.V_G6PD_uM_min, V_PGLS = row.V_PGLS_uM_min,
       V_PGD = row.V_PGD_uM_min, V_RPI = row.V_RPI_uM_min,
       V_RPE = row.V_RPE_uM_min, V_TKT_Rxn1 = row.V_TKT_Rxn1_uM_min,
       V_TKT_Rxn2 = row.V_TKT_Rxn2_uM_min, V_TA = row.V_TA_uM_min,
       V_R5Pase = row.V_R5Pase_uM_min)

"""
    row_to_pools(row; order = _PANEL_METAB_ORDER) -> Vector{Float64}

Pull steady-state pool concentrations (µM) for the panel metabolites out of one grid-CSV
`DataFrameRow`, reading the `<sym>_uM` columns in `order`. NaN pool columns (masked rows)
pass through as NaN.
"""
row_to_pools(row; order = _PANEL_METAB_ORDER) =
    Float64[getproperty(row, Symbol(string(sym), "_uM")) for sym in order]

"""
    nearest_cell(xs, ys, xc, yc) -> (i_r5p, i_nadph)

Snap a click at data coords `(xc, yc)` to the nearest grid cell on the log-spaced demand
axes `xs` (sorted unique `r5p_phi`) and `ys` (sorted unique `nadph_phi`). Distance is
measured in log10 space (the heatmap axes are log10-scaled).
"""
nearest_cell(xs, ys, xc, yc) =
    (argmin(abs.(log10.(xs) .- log10(xc))), argmin(abs.(log10.(ys) .- log10(yc))))

"""
    log_cell_edges(centers) -> Vector  (length = length(centers) + 1)

Cell boundaries for a `heatmap!` whose axis is log10-scaled. Given the log-evenly-spaced cell
centers, return the n+1 edges sitting midway between centers IN LOG SPACE, outer two
extrapolated by half a log-step. Keeps cells equal-sized on the log axis.
"""
function log_cell_edges(centers)
    lc = log10.(centers)
    d  = length(lc) > 1 ? (lc[end] - lc[1]) / (length(lc) - 1) : 1.0
    return 10 .^ vcat(lc .- d / 2, lc[end] + d / 2)
end

"""
    load_grid(path) -> DataFrame

Read a grid CSV and normalize the `retcode` column from String to Symbol (CSV round-trips it
as "Terminated"; the renderers mask on the `:Terminated` Symbol).
"""
function load_grid(path)
    df = CSV.read(path, DataFrame)
    df.retcode = Symbol.(df.retcode)
    return df
end
