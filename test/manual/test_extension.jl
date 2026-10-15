# Manual extension check. PPP.jl is unregistered, so CI cannot load the extension; run this by
# hand from the PPP_Experiments environment, which `develop`s this package and PPP.jl by path:
#   julia --project=/home/james/projects/PPP_Experiments test/manual/test_extension.jl
using Test, CSV, DataFrames
using GlycolysisPentosePathwayPhaseMap, PentosePhosphatePathway
const G   = GlycolysisPentosePathwayPhaseMap
const PPP = PentosePhosphatePathway

ext = Base.get_extension(G, :PentosePhosphatePathwayExt)

@testset "recompute extension" begin
    @test ext !== nothing

    # initial state = PPP.jl's canonical combined init, except F26BP at the phase map's
    # physiological 5 µM (cytosolic water; PPP.jl's own default is 0)
    s = ext._demand_setup(0.10, 2, 2, 1e-3, 1.0)
    @test keys(s.init0) == keys(PPP.glycolysis_ppp_init_conc)
    @test s.init0.F26BP == G._PHYSIOLOGICAL_DEMAND.f26bp_uM * 1e-6
    @test s.init0.F26BP ≈ 5e-6 rtol = 1e-12      # 5.0 * 1e-6 is one ulp below the literal 5e-6
    same = [k for k in keys(s.init0) if k != :F26BP]
    @test all(getproperty(s.init0, k) == getproperty(PPP.glycolysis_ppp_init_conc, k) for k in same)
    @test s.init0 !== PPP.glycolysis_ppp_init_conc          # a copy, never the package global
    @test PPP.glycolysis_ppp_init_conc.F26BP == 0.0         # the package global is untouched
    s0 = ext._demand_setup(0.10, 2, 2, 1e-3, 1.0; f26bp_uM = 0.0)
    @test collect(s0.init0) == collect(PPP.glycolysis_ppp_init_conc)   # 0 reproduces PPP.jl

    # 2 x 2 smoke regenerate into scratch paths, with provenance
    dir = mktempdir()
    out = regenerate_grid(; atp_out = joinpath(dir, "atp.csv"),
                            demand_out = joinpath(dir, "demand.csv"),
                            n_nadph = 2, n_r5p = 2, atp_fracs = [0.10])
    @test out.provenance_out == joinpath(dir, "provenance.toml")
    atp = G.load_grid(out.atp_out)
    @test nrow(atp) == 4
    @test count(==(:Terminated), atp.retcode) >= 1
    prov = G.read_provenance(out.provenance_out)["grids"]
    gs = G.git_state(pkgdir(PPP))
    @test prov["ppp_commit"] == gs.commit
    @test prov["ppp_branch"] == gs.branch
    @test prov["ppp_dirty"] == gs.dirty
    @test prov["atp_grid_rows"] == 4
    @test prov["demand_grid_rows"] == 4
    @test prov["atp_grid_nonconverged"] == count(!=(:Terminated), atp.retcode)
    @test prov["atp_fracs"] == [0.10]
    @test prov["atpase_frac"] == G._PHYSIOLOGICAL_DEMAND.atpase_frac   # default 2-D background
    @test prov["f26bp_uM"] == G._PHYSIOLOGICAL_DEMAND.f26bp_uM
    for s in (:ADP, :AMP, :NAD, :NADH)
        @test "$(s)_uM" in names(atp)
    end
    @test !haskey(G.read_provenance(out.provenance_out), "literature")
end

@testset "cold start, flux units, physiological defaults" begin
    # PhaseMap's physiological demand constants are PPP.jl's default sink settings, as fractions
    # of the same supply ceilings the grids use
    s1 = ext._demand_setup(0.25, 2, 2, 1e-3, 1.0)
    @test G._PHYSIOLOGICAL_DEMAND.nadph_phi ≈ PPP.PPP_params.NADPHox_Vmax / s1.S_NADPH rtol = 1e-6
    @test G._PHYSIOLOGICAL_DEMAND.r5p_phi ≈ PPP.PPP_params.R5Pase_Vmax / s1.S_R5P rtol = 1e-6

    # cached fluxes are µM/min (model rates are M/min): the R5P sink at PPP.jl's default φ
    # delivers ≈ its 38 µM/min anchor (Fan 2014), not 60× that
    s = ext._demand_setup(0.25, 1, 1, 0.059, 0.059)
    r = only(ext._demand_column(s, 1, 0.1068))
    @test r.retcode == :Terminated
    @test 30 < r.V_R5Pase_uM_min < 40

    # cold start: a cell's steady state does not depend on the cells solved before it in its
    # column (a warm-started sweep dragged this whole column into ATP collapse at 0.5·HK1_cap)
    s12 = ext._demand_setup(0.25, 12, 12, 1e-3, 1.0)
    col = ext._demand_column(s12, 8, s12.r5p_phis[8])
    @test col[7].retcode == :Terminated && col[7].ATP_uM > 1000
    alone = ext._demand_setup(0.25, 1, 1, s12.nadph_phis[7], s12.nadph_phis[7])
    @test only(ext._demand_column(alone, 8, s12.r5p_phis[8])).ATP_uM ≈ col[7].ATP_uM rtol = 1e-6

    # the four cofactor pools are cached, finite and positive at the physiological cell
    # (AMP included: a vanished AMP is what the panel's below-floor arrow is for)
    for k in (:ADP_uM, :AMP_uM, :NAD_uM, :NADH_uM)
        v = getproperty(r, k)
        @test isfinite(v) && v > 0
    end
    # F26BP sets the adenylate split (spec D7): 5 µM keeps ADP in the free-ADP data range
    @test 10 < r.ADP_uM < 100
end

