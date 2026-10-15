##########################################################################################
#        PentosePhosphatePathwayExt.jl — live grid recompute (kinetic-model path)        #
##########################################################################################
#
# Activated only when PentosePhosphatePathway (and the solver stack it pulls in) is loaded.
# Provides the real recompute implementation (`_regenerate_grid_impl`) and registers it into the
# core module's `_REGEN_HOOK` Ref via __init__ (NOT by overriding the parent's regenerate_grid
# method, which would break precompilation): it solves the integrated glycolysis+PPP ODE model
# over the demand grids and rewrites the cached CSVs the viewer reads.
#
# Sweep logic was ported from the origin PPP_flux_simulator_core.jl, with two corrections: every
# cell is solved cold from PPP.jl's canonical initial state (a warm-started sweep carried an
# ATP-collapsed state down whole columns), and fluxes are cached as µM/min from the model's M/min
# (×1e6; the port's extra ×60 assumed per-second rates). The pure diagnostics
# pentose_cycle_index/classify_mode live in the core module's grid_io.jl.

module PentosePhosphatePathwayExt

using GlycolysisPentosePathwayPhaseMap
using PentosePhosphatePathway, OrdinaryDiffEq, DiffEqCallbacks, LabelledArrays
using DataFrames
import CSV
import CellMetabolism, CellMetabolismBase

const PPP  = PentosePhosphatePathway
const GPPM = GlycolysisPentosePathwayPhaseMap

# Pools / fluxes cached per cell (mirrors the viewer's expected CSV columns).
const _DEMAND_POOLS = (:R5P, :S7P, :NADPH, :NADP, :Ru5P, :X5P, :E4P, :F6P, :GAP, :G6P,
                       :Glucose, :F16BP, :DHAP, :PGLn, :PGA, :ATP, :ADP, :AMP, :NAD, :NADH)
const _DEMAND_FLUXK = (:V_G6PD, :V_PGLS, :V_PGD, :V_RPI, :V_RPE, :V_TKT_Rxn1, :V_TKT_Rxn2,
                       :V_TA, :V_R5Pase, :V_NADPHox,
                       :V_HK1, :V_GPI, :V_PFKP, :V_ALDO, :V_TPI, :V_GAPDH)

# Solve the integrated model to steady state with the repo-standard solver settings.
function steady_state(params; u0)
    prob = make_ODEProblem(glycolysis_ppp_pathway, u0, (0.0, 1e10), params)
    solve(prob, Rodas5P();
          callback = TerminateSteadyState(1e-13, 1e-8),
          abstol = 1e-14, reltol = 1e-9)
end

# Shared per-sweep state: base params/init, the NADPH/R5P supply ceilings, the demand-φ axes.
function _demand_setup(atpase_frac, n_nadph, n_r5p, phi_lo, phi_hi;
                       f26bp_uM = GPPM._PHYSIOLOGICAL_DEMAND.f26bp_uM)
    params0 = PPP.merge_LArrays(PPP.PPP_params, CellMetabolism.glycolysis_params)
    # PPP.jl's canonical combined init: PPP_init_conc, then PPP_glycolysis_init_overrides
    # (FIX-15/16/22), then CellMetabolism's. A hand-merge skips the overrides. Copied: the
    # export is a package-level array that the sweep must never write into.
    init0   = copy(PPP.glycolysis_ppp_init_conc)
    # F26BP is cytosolic water (no whole_cell_to_cytosolic), as in PPP_Experiments'
    # Exp_basic_panel_v2.jl. It has no reactions, so it holds this value at steady state.
    init0.F26BP = f26bp_uM * 1e-6
    RPI_cap = params0.RPI_Vmax * params0.RPI_Conc
    S_NADPH = 6 * RPI_cap                  # 2·V_ox,max, V_ox,max = 3·RPI_cap (Ru5P-disposal limit)
    S_R5P   = 1.5 * RPI_cap                # RPI fwd + TKT_Rxn1-reverse bypass; 2/3 of R5P via RPI
    params0.ATPase_Vmax = atpase_frac * (params0.HK1_Vmax * params0.HK1_Conc * 2)  # φ·S_ATP
    nadph_phis = 10 .^ range(log10(phi_lo), log10(phi_hi), length = n_nadph)
    r5p_phis   = 10 .^ range(log10(phi_lo), log10(phi_hi), length = n_r5p)
    return (; params0, init0, S_NADPH, S_R5P, nadph_phis, r5p_phis)
end

# Solve one R5P-demand column down the NADPH-demand axis. Every cell starts cold from init0, so
# its steady state never depends on sweep order (the model has an absorbing ATP-collapsed state
# that a warm start would carry from one cell into the next).
function _demand_column(setup, ir, rphi)
    (; params0, init0, S_NADPH, S_R5P, nadph_phis) = setup
    rows = NamedTuple[]
    for (in_, nphi) in enumerate(nadph_phis)
        p = copy(params0)
        p.NADPHox_Vmax = nphi * S_NADPH
        p.R5Pase_Vmax  = rphi * S_R5P

        sol = try
            steady_state(p; u0 = copy(init0))
        catch
            nothing
        end
        converged = sol !== nothing && Symbol(sol.retcode) == :Terminated
        u = converged ? sol.u[end] : (sol === nothing ? nothing : sol.u[end])

        if u !== nothing
            f  = PPP.calc_all_fluxes(u, p)
            ci = GPPM.pentose_cycle_index(f)
            md = GPPM.classify_mode(f)
            row = merge(
                (; i_nadph = in_, i_r5p = ir,
                   nadph_phi = nphi, r5p_phi = rphi,
                   NADPHox_Vmax = p.NADPHox_Vmax, R5Pase_Vmax = p.R5Pase_Vmax,
                   retcode = sol === nothing ? :SolverError : Symbol(sol.retcode),
                   cycle_index = ci, mode = md),
                NamedTuple{Tuple(Symbol.(string.(_DEMAND_FLUXK), "_uM_min"))}(
                    Tuple(getfield(f, k) * 1e6 for k in _DEMAND_FLUXK)),
                NamedTuple{Tuple(Symbol.(string.(_DEMAND_POOLS), "_uM"))}(
                    Tuple(getproperty(u, k) * 1e6 for k in _DEMAND_POOLS)),
            )
            push!(rows, row)
        else
            push!(rows, (; i_nadph = in_, i_r5p = ir,
                nadph_phi = nphi, r5p_phi = rphi,
                NADPHox_Vmax = p.NADPHox_Vmax, R5Pase_Vmax = p.R5Pase_Vmax,
                retcode = :SolverError, cycle_index = NaN, mode = :undetermined,
                (Symbol(string(k), "_uM_min") => NaN for k in _DEMAND_FLUXK)...,
                (Symbol(string(k), "_uM") => NaN for k in _DEMAND_POOLS)...))
        end
    end
    return rows
end

# 2-D NADPH × R5P demand grid at a fixed ATP background.
function run_demand_grid(; n_nadph = 12, n_r5p = 12,
                           atpase_frac = GPPM._PHYSIOLOGICAL_DEMAND.atpase_frac,
                           phi_lo = 1e-3, phi_hi = 1.0,
                           f26bp_uM = GPPM._PHYSIOLOGICAL_DEMAND.f26bp_uM)
    setup = _demand_setup(atpase_frac, n_nadph, n_r5p, phi_lo, phi_hi; f26bp_uM)
    rows = NamedTuple[]
    for (ir, rphi) in enumerate(setup.r5p_phis)
        append!(rows, _demand_column(setup, ir, rphi))
    end
    return DataFrame(rows)
end

# ATP-demand levels of the 3-D grid: 11 log-spaced fractions of the 2·HK1_cap ATP ceiling from
# 0.01 up to the physiological operating point's 0.25 (= 0.5·HK1_cap).
_default_atp_fracs() =
    10 .^ range(log10(0.01), log10(GPPM._PHYSIOLOGICAL_DEMAND.atpase_frac), length = 11)

# 3-D ATP × NADPH × R5P grid: the demand grid repeated across log-spaced ATP-demand levels,
# each block tagged with its `atpase_frac`. Serial (run-once, offline) — minutes for a full grid.
function run_atp_grid(; n_nadph = 12, n_r5p = 12, atp_fracs = _default_atp_fracs(),
                        phi_lo = 1e-3, phi_hi = 1.0,
                        f26bp_uM = GPPM._PHYSIOLOGICAL_DEMAND.f26bp_uM)
    frames = DataFrame[]
    for f in atp_fracs
        s = _demand_setup(f, n_nadph, n_r5p, phi_lo, phi_hi; f26bp_uM)
        for ir in 1:n_r5p
            sub = DataFrame(_demand_column(s, ir, s.r5p_phis[ir]))
            sub.atpase_frac .= f
            push!(frames, sub)
        end
    end
    return vcat(frames...)
end

# Real implementation behind GPPM.regenerate_grid. Registered into the core module's hook Ref by
# __init__ (below) at load time — never defined as a method on GPPM.regenerate_grid, which would
# be a precompile-breaking method overwrite. Docstring lives on the exported parent function.
function _regenerate_grid_impl(; atp_out = GPPM.default_atp_grid(),
                                 demand_out = GPPM.default_demand_grid(),
                                 provenance_out = joinpath(dirname(atp_out), "provenance.toml"),
                                 n_nadph = 12, n_r5p = 12,
                                 atp_fracs = _default_atp_fracs(),
                                 atpase_frac = GPPM._PHYSIOLOGICAL_DEMAND.atpase_frac,
                                 phi_lo = 1e-3, phi_hi = 1.0,
                                 f26bp_uM = GPPM._PHYSIOLOGICAL_DEMAND.f26bp_uM)
    mkpath(dirname(atp_out)); mkpath(dirname(demand_out))
    atp = run_atp_grid(; n_nadph, n_r5p, atp_fracs, phi_lo, phi_hi, f26bp_uM)
    CSV.write(atp_out, atp)
    demand = run_demand_grid(; n_nadph, n_r5p, atpase_frac, phi_lo, phi_hi, f26bp_uM)
    CSV.write(demand_out, demand)
    GPPM.update_provenance!(provenance_out, "grids",
        _grid_provenance(atp, demand; n_nadph, n_r5p, atp_fracs, atpase_frac, phi_lo, phi_hi,
                         f26bp_uM))
    @info "Regenerated grids" atp_out demand_out provenance_out atp_rows=nrow(atp) demand_rows=nrow(demand)
    return (; atp_out, demand_out, provenance_out)
end

# The [grids] provenance table: which model build produced the grids, with what axes, and how
# many cells failed to reach steady state. A git checkout of PPP.jl records commit/branch/dirty;
# a registry install records its version instead.
function _grid_provenance(atp, demand; n_nadph, n_r5p, atp_fracs, atpase_frac, phi_lo, phi_hi,
                          f26bp_uM)
    gs = GPPM.git_state(pkgdir(PPP))
    prov = Dict{String,Any}(
        "ppp_commit"               => gs === nothing ? "v" * string(pkgversion(PPP)) : gs.commit,
        "cellmetabolism_version"   => string(pkgversion(CellMetabolism)),
        "generated"                => Libc.strftime("%Y-%m-%d", time()),
        "n_nadph" => n_nadph, "n_r5p" => n_r5p,
        "atp_fracs" => collect(Float64, atp_fracs), "atpase_frac" => Float64(atpase_frac),
        "phi_lo" => Float64(phi_lo), "phi_hi" => Float64(phi_hi),
        "f26bp_uM" => Float64(f26bp_uM),
        "atp_grid_rows"            => nrow(atp),
        "demand_grid_rows"         => nrow(demand),
        "atp_grid_nonconverged"    => count(!=(:Terminated), atp.retcode),
        "demand_grid_nonconverged" => count(!=(:Terminated), demand.retcode),
    )
    if gs !== nothing
        prov["ppp_branch"] = gs.branch
        prov["ppp_dirty"]  = gs.dirty
    end
    return prov
end

__init__() = (GPPM._REGEN_HOOK[] = _regenerate_grid_impl)

end # module
