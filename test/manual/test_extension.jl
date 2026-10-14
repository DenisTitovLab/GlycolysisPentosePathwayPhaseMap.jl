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

    # initial state = PPP.jl's canonical combined init (with PPP_glycolysis_init_overrides)
    s = ext._demand_setup(0.10, 2, 2, 1e-3, 1.0)
    @test keys(s.init0) == keys(PPP.glycolysis_ppp_init_conc)
    @test collect(s.init0) == collect(PPP.glycolysis_ppp_init_conc)
    @test s.init0 !== PPP.glycolysis_ppp_init_conc          # a copy, never the package global

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
end

