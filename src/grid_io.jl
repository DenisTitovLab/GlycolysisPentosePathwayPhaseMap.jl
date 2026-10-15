##########################################################################################
#                  grid_io.jl — pure flux diagnostics + grid-CSV helpers                 #
##########################################################################################

# Metabolite panel, split as in PPP_Experiments' basic panel. Left: species whose free and
# total pools match, compared with extraction metabolomics. Right: cofactors whose pool is
# mostly protein-bound, compared with free-pool data only. The grouping is PPP_Experiments'
# scripts/helper_files/metab_data.jl (FREE_EQ_TOTAL_METABS / FREE_POOL_METABS); its literature
# export (check_panel_grouping) refuses to run if these tuples disagree with it.
const _FREE_EQ_TOTAL_PANEL = (:Glucose, :G6P, :F6P, :F16BP, :DHAP, :GAP, :ATP,
                              :PGLn, :PGA, :Ru5P, :R5P, :X5P, :S7P, :E4P)
const _FREE_POOL_PANEL     = (:ADP, :AMP, :NAD, :NADH, :NADP, :NADPH)
# Display order = the two panels concatenated. Single source of truth for row_to_pools,
# load_metab_lit, the metabolite panel and the PPP_Experiments export.
const _PANEL_METAB_ORDER   = (_FREE_EQ_TOTAL_PANEL..., _FREE_POOL_PANEL...)

# Default operating point: PPP.jl's default sink settings (src/model_parameters.jl:
# NADPHox_Vmax = 0.059·(6·RPI_cap), R5Pase_Vmax = 0.1068·(1.5·RPI_cap)) as fractions φ of the
# same supply ceilings the grids use, at an ATPase demand of 0.5·HK1_cap (= 0.25 of the
# 2·HK1_cap ATP ceiling; PPP.jl's own default is 0.2·HK1_cap). F26BP is 5 µM (cytosolic water),
# the physiological operating point of PPP_Experiments
# notes/2026-09-23_metabolite-dynamics-physiological-operating-point.md and the value its basic
# panels use; PPP.jl's own default is 0, at which upper glycolysis runs near its limit and ADP
# builds up ~20-fold. The 2-D grid is solved at this ATP level, and phasemap()/explore() open
# on the grid cell nearest it. The manual extension test pins the two φ against PPP.jl.
const _PHYSIOLOGICAL_DEMAND = (nadph_phi = 0.059, r5p_phi = 0.1068, atpase_frac = 0.25,
                               f26bp_uM = 5.0)

"""
    pentose_cycle_index(fluxes; tol = 1e-12) -> Float64

Fraction of oxidative-PPP-generated pentose carbon recycled rather than exported as
R5P:  (V_ox - V_R5Pase) / V_ox  with  V_ox = fluxes.V_G6PD.
≈1 full pentose cycle (max NADPH/glucose); ≈0 linear; <0 net non-oxidative R5P
production (non-ox runs reverse). Returns NaN when V_ox ≤ tol (ill-defined) or when the
R5P sink runs backward (V_R5Pase < 0, which only happens once ATP has collapsed): the index
would exceed 1, which no recycled fraction can.
"""
function pentose_cycle_index(fluxes; tol = 1e-12)
    V_ox = fluxes.V_G6PD
    (V_ox <= tol || fluxes.V_R5Pase < 0) && return NaN
    return (V_ox - fluxes.V_R5Pase) / V_ox
end

"""
    classify_mode(fluxes; cycle_thresh = 0.5) -> Symbol

:cycle (index ≥ cycle_thresh), :reverse (index ≤ 0), :linear (in between), or
:undetermined (index is NaN: oxPPP flux ≈ 0, or the R5P sink running backward).
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
    physiological_cell(df) -> (i_nadph, i_r5p)

The grid cell nearest (in log space) to the physiological demand `_PHYSIOLOGICAL_DEMAND`:
PPP.jl's default NADPH and R5P sink settings. Keyed like the explorer's `selected`.
"""
function physiological_cell(df)
    xs = sort(unique(df.r5p_phi)); ys = sort(unique(df.nadph_phi))
    ir, in_ = nearest_cell(xs, ys, _PHYSIOLOGICAL_DEMAND.r5p_phi, _PHYSIOLOGICAL_DEMAND.nadph_phi)
    return (in_, ir)
end

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

# Columns every renderer reads (phasemap and explore). Flux columns come from _PPP_FLUX_KEYS
# (flux_map.jl, included after this file; resolved at call time).
const _GRID_CORE_COLUMNS = ("retcode", "cycle_index", "mode", "i_nadph", "i_r5p",
                            "nadph_phi", "r5p_phi")

"Throw an `ArgumentError` naming every entry of `cols` missing from `df`, with reason `why`."
function _require_columns(df, cols, why)
    absent = [c for c in cols if !(c in names(df))]
    isempty(absent) || throw(ArgumentError("grid is missing column(s) $(join(absent, ", ")): $why"))
    return df
end

"""
    load_grid(path) -> DataFrame

Read a grid CSV, check it has the columns every renderer needs (`_GRID_CORE_COLUMNS` and a
`<key>_uM_min` flux column per `_PPP_FLUX_KEYS`), and normalize `retcode` from String to Symbol
(CSV round-trips it as "Terminated"; the renderers mask on the `:Terminated` Symbol). Pool
columns are checked by `build_explorer`, the only reader.
"""
function load_grid(path)
    df = CSV.read(path, DataFrame)
    _require_columns(df, [_GRID_CORE_COLUMNS..., ("$(k)_uM_min" for k in _PPP_FLUX_KEYS)...],
                     "not a phase-map grid CSV")
    df.retcode = Symbol.(df.retcode)
    return df
end
