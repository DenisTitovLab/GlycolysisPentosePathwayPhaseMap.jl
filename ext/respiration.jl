##########################################################################################
#              respiration.jl — ATPSynthase respiration pseudo-enzyme (what-if)          #
##########################################################################################
#
# Respiration modelled as one ATP-producing reaction, ADP + Pi -> ATP, following Kober et al.
# 2026 (bioRxiv 10.64898/2026.01.31.703051, Methods "Glycolysis Model Code"; code in
# CellMetabolism.jl examples/2026_F26BP_Figures/). It has no carbon or NADH cost: respiration
# competes with glycolysis only for ADP and Pi.
#
# Kober's construction holds the total ATP-supply capacity fixed and splits it by a respiration
# fraction r: ATPSynthase_Vmax = r · C and every glycolytic enzyme concentration × (1 − r), with
# C = 0.3 · 2·HK1_cap (glycolysis delivers at most ~30% of its 2·HK1_cap net-yield bound).
# r = 0 is the unmodified model. PPP enzymes, sinks, kinetic constants and inits are untouched.

# Rate-law constants of Kober's Methods (M; Keq in M⁻¹). Model choices with no literature
# source: Km_ATP = 1000 M switches off product inhibition, and Km_ADP = Km_Pi = 10 µM.
const _ATPSYNTHASE_KM_ATP       = 1000.0
const _ATPSYNTHASE_KM_ADP       = 1e-5
const _ATPSYNTHASE_KM_PHOSPHATE = 1e-5
const _ATPSYNTHASE_KEQ          = 1e6

# C = _RESP_CAPACITY_FACTOR · 2·HK1_cap, from the unscaled (r = 0) HK1.
const _RESP_CAPACITY_FACTOR = 0.3

# Kober's scaled set, restricted to the enzymes with a concentration parameter in the coupled
# model (AK and NDPK have none).
const _GLYCOLYTIC_CONC = (:GLUT1_Conc, :HK1_Conc, :GPI_Conc, :PFKP_Conc, :ALDO_Conc, :TPI_Conc,
                          :GAPDH_Conc, :PGK_Conc, :PGAM_Conc, :ENO_Conc, :PKM2_Conc, :LDH_Conc,
                          :MCT_Conc)

const _ATPSYNTHASE_SPEC = (:ATPSynthase, (:ADP, :Phosphate), (:ATP,))

# A pathway with one more reaction, built from the type parameters of an existing one so it
# tracks PPP.jl's reaction list instead of copying it.
_append_enzyme(::CellMetabolismBase.MetabolicPathway{C,E}, spec) where {C,E} =
    CellMetabolismBase.MetabolicPathway(C, (E..., spec))

const _RESP_PATHWAY = _append_enzyme(PPP.glycolysis_ppp_pathway, _ATPSYNTHASE_SPEC)

# Reversible Michaelis–Menten with competitive product inhibition (Kober Methods): the
# CellMetabolism ATPase law with the direction reversed.
@inline function CellMetabolismBase.rate(
    ::CellMetabolismBase.Enzyme{:ATPSynthase,(:ADP, :Phosphate),(:ATP,)},
    metabs,
    params,
)
    return (params.ATPSynthase_Vmax / (_ATPSYNTHASE_KM_ADP * _ATPSYNTHASE_KM_PHOSPHATE)) *
           (1 / (1 +
                 metabs.ATP / _ATPSYNTHASE_KM_ATP +
                 metabs.ADP / _ATPSYNTHASE_KM_ADP +
                 metabs.Phosphate / _ATPSYNTHASE_KM_PHOSPHATE +
                 (metabs.ADP / _ATPSYNTHASE_KM_ADP) *
                 (metabs.Phosphate / _ATPSYNTHASE_KM_PHOSPHATE))) *
           (metabs.ADP * metabs.Phosphate - metabs.ATP / params.ATPSynthase_Keq)
end

"""
    _with_respiration(params0, resp_frac) -> LArray

`params0` with a respiration fraction `resp_frac` of the fixed total ATP-supply capacity moved
from glycolysis to ATPSynthase: the `_GLYCOLYTIC_CONC` concentrations × (1 − resp_frac), and
`ATPSynthase_Vmax = resp_frac · _RESP_CAPACITY_FACTOR · 2·HK1_cap` (unscaled HK1_cap).
"""
function _with_respiration(params0, resp_frac)
    0 <= resp_frac < 1 || throw(ArgumentError("resp_frac must be in [0, 1), got $resp_frac"))
    HK1_cap = params0.HK1_Vmax * params0.HK1_Conc
    p = PPP.merge_LArrays(copy(params0), LVector(
        ATPSynthase_Vmax = resp_frac * _RESP_CAPACITY_FACTOR * 2 * HK1_cap,
        ATPSynthase_Keq  = _ATPSYNTHASE_KEQ))
    for k in _GLYCOLYTIC_CONC
        setproperty!(p, k, getproperty(p, k) * (1 - resp_frac))
    end
    return p
end

# ATPSynthase flux (M/min). PPP.calc_all_fluxes enumerates only glycolysis_ppp_pathway.
_resp_flux(u, p) = CellMetabolismBase.rate(CellMetabolismBase.Enzyme(_ATPSYNTHASE_SPEC...), u, p)
