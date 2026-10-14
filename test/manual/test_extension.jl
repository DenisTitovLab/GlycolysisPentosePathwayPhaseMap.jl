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
    @test !haskey(G.read_provenance(out.provenance_out), "literature")
end
