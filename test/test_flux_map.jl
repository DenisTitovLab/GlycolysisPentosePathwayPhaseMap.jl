using Test, CairoMakie, GlycolysisPentosePathwayPhaseMap
const G = GlycolysisPentosePathwayPhaseMap
CairoMakie.activate!()

@testset "draw_ppp_flux_map!" begin
    @test haskey(G._PPP_NODES, :GAP_p)
    @test :V_TKT_Rxn2 in G._PPP_FLUX_KEYS
    base = (; V_HK1 = 9e-6, V_GPI = 8e-6, V_PFKP = 8e-6, V_ALDO = 8e-6, V_TPI = 4e-6,
              V_GAPDH = 8e-6, V_G6PD = 6e-7, V_PGL = 6e-7, V_PGD = 6e-7, V_RPI = 2e-7,
              V_RPE = 4e-7, V_TKT_Rxn1 = 2e-7, V_TKT_Rxn2 = 2e-7, V_TA = 2e-7, V_R5Pase = 1e-7)
    fig = Figure(size = (500, 460)); ax = Axis(fig[1, 1])
    @test_nowarn G.draw_ppp_flux_map!(ax, base; title = "fwd")
    out = joinpath(tempdir(), "ppp_flux_map_smoke.png"); save(out, fig); @test isfile(out)
    rev = merge(base, (; V_TKT_Rxn1 = -2e-7, V_TA = -2e-7, V_RPE = -4e-7))
    fig2 = Figure(size = (500, 460)); ax2 = Axis(fig2[1, 1])
    @test_nowarn G.draw_ppp_flux_map!(ax2, rev; title = "rev")
end
