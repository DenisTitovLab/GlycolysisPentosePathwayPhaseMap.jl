##########################################################################################
#              respiration.jl — ATPSynthase respiration pseudo-enzyme (what-if)          #
##########################################################################################
#
# Respiration modelled as one carbon-linked ATP-producing reaction,
#     Pyruvate + NADH + 15 ADP + 15 Pi -> NAD + 3 CO2 + 15 ATP,
# full pyruvate oxidation (PDH + TCA + OxPhos, 12.5 ATP) plus the cytosolic NADH that glycolysis
# made for it, carried in by the malate-aspartate shuttle (2.5 ATP); with glycolysis's 2 ATP
# that is 32 ATP per glucose. O2 and H2O are implicit; CO2 is the model's fixed boundary species.
# Consuming the NADH is what lets the reaction carry flux at steady state: the model conserves
# BPG + 3PG + 2PG + PEP + Pyr − NADH, so a pyruvate-only sink would leave LDH as the sole NADH
# reoxidizer and force every pyruvate through it.
#
# The capacity construction follows Kober et al. 2026 (bioRxiv 10.64898/2026.01.31.703051,
# Methods "Glycolysis Model Code"; code in CellMetabolism.jl examples/2026_F26BP_Figures/), who
# modelled respiration as the carbon-free ADP + Pi -> ATP: the total ATP-supply capacity is held
# fixed and split by a respiration fraction r. Respiration gets ATP capacity r · C
# (ATPSynthase_Vmax = r · C / 15, in pyruvate units) and every glycolytic enzyme concentration is
# scaled by (1 − r), with C = 0.3 · 2·HK1_cap (glycolysis delivers at most ~30% of its 2·HK1_cap
# net-yield bound). r = 0 is the unmodified model. PPP enzymes, sinks, kinetic constants and
# inits are untouched.

# ATP made per turnover (= per pyruvate oxidized).
const _RESP_ATP_YIELD = 15

# Rate-law constants (M). The ADP and Pi saturation term is Kober's (model choices with no
# literature source; Km_ATP = 1000 M switches off product inhibition). K_PYRUVATE is LDH's Km
# for pyruvate (CellMetabolism glycolysis_params; the mitochondrial pyruvate carrier's is
# similar). K_NADH has no source: it is set ~10x below the model's free cytosolic NADH
# (~0.22 µM at the physiological cell) so respiration is ADP-controlled at normal redox and
# shuts off only when NADH or pyruvate drain. LDH's own Km_NADH (20 µM) would starve it.
const _ATPSYNTHASE_KM_ATP       = 1000.0
const _ATPSYNTHASE_KM_ADP       = 1e-5
const _ATPSYNTHASE_KM_PHOSPHATE = 1e-5
const _ATPSYNTHASE_K_PYRUVATE   = 140e-6
const _ATPSYNTHASE_K_NADH       = 2e-8

# C = _RESP_CAPACITY_FACTOR · 2·HK1_cap, from the unscaled (r = 0) HK1.
const _RESP_CAPACITY_FACTOR = 0.3

# Kober's scaled set, restricted to the enzymes with a concentration parameter in the coupled
# model (AK and NDPK have none).
const _GLYCOLYTIC_CONC = (:GLUT1_Conc, :HK1_Conc, :GPI_Conc, :PFKP_Conc, :ALDO_Conc, :TPI_Conc,
                          :GAPDH_Conc, :PGK_Conc, :PGAM_Conc, :ENO_Conc, :PKM2_Conc, :LDH_Conc,
                          :MCT_Conc)

const _ATPSYNTHASE_SPEC = (:ATPSynthase,
    (:Pyruvate, :NADH, ntuple(_ -> :ADP, _RESP_ATP_YIELD)..., ntuple(_ -> :Phosphate, _RESP_ATP_YIELD)...),
    (:NAD, :CO2, :CO2, :CO2, ntuple(_ -> :ATP, _RESP_ATP_YIELD)...))

# A pathway with one more reaction, built from the type parameters of an existing one so it
# tracks PPP.jl's reaction list instead of copying it.
_append_enzyme(::CellMetabolismBase.MetabolicPathway{C,E}, spec) where {C,E} =
    CellMetabolismBase.MetabolicPathway(C, (E..., spec))

const _RESP_PATHWAY = _append_enzyme(PPP.glycolysis_ppp_pathway, _ATPSYNTHASE_SPEC)

# Irreversible respiratory control: Kober's ADP·Pi saturation term (same denominator, no
# reverse term) times hyperbolic pyruvate and NADH supply terms. Rate is per pyruvate.
@inline function CellMetabolismBase.rate(::CellMetabolismBase.Enzyme{:ATPSynthase}, metabs, params)
    adp = metabs.ADP / _ATPSYNTHASE_KM_ADP
    pi_ = metabs.Phosphate / _ATPSYNTHASE_KM_PHOSPHATE
    return params.ATPSynthase_Vmax *
           (adp * pi_ / (1 + metabs.ATP / _ATPSYNTHASE_KM_ATP + adp + pi_ + adp * pi_)) *
           (metabs.Pyruvate / (_ATPSYNTHASE_K_PYRUVATE + metabs.Pyruvate)) *
           (metabs.NADH / (_ATPSYNTHASE_K_NADH + metabs.NADH))
end

"""
    _with_respiration(params0, resp_frac) -> LArray

`params0` with a respiration fraction `resp_frac` of the fixed total ATP-supply capacity moved
from glycolysis to ATPSynthase: the `_GLYCOLYTIC_CONC` concentrations × (1 − resp_frac), and
`ATPSynthase_Vmax = resp_frac · _RESP_CAPACITY_FACTOR · 2·HK1_cap / _RESP_ATP_YIELD` (unscaled
HK1_cap; pyruvate units, so the ATP capacity is resp_frac · C).
"""
function _with_respiration(params0, resp_frac)
    0 <= resp_frac < 1 || throw(ArgumentError("resp_frac must be in [0, 1), got $resp_frac"))
    HK1_cap = params0.HK1_Vmax * params0.HK1_Conc
    p = PPP.merge_LArrays(copy(params0), LVector(
        ATPSynthase_Vmax = resp_frac * _RESP_CAPACITY_FACTOR * 2 * HK1_cap / _RESP_ATP_YIELD))
    for k in _GLYCOLYTIC_CONC
        setproperty!(p, k, getproperty(p, k) * (1 - resp_frac))
    end
    return p
end

# ATPSynthase flux (M/min of pyruvate; ATP is _RESP_ATP_YIELD × this). PPP.calc_all_fluxes enumerates only glycolysis_ppp_pathway.
_resp_flux(u, p) = CellMetabolismBase.rate(CellMetabolismBase.Enzyme(_ATPSYNTHASE_SPEC...), u, p)
