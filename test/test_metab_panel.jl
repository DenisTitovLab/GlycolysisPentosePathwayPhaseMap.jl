using Test, CairoMakie, DataFrames, GlycolysisPentosePathwayPhaseMap
const G = GlycolysisPentosePathwayPhaseMap
CairoMakie.activate!()

@testset "load_metab_lit" begin
    lit = G.load_metab_lit(G.metab_lit_path())
    @test Set(keys(lit)) == Set(G._PANEL_METAB_ORDER)
    @test eltype(lit[:G6P]) == Float64
    @test !isempty(lit[:G6P]) && all(>(0), lit[:G6P])
    @test lit[:PGLn] isa Vector{Float64}
end

@testset "draw_metab_panel! offscreen" begin
    lit = Dict(s => Float64[] for s in G._PANEL_METAB_ORDER)
    lit[:G6P] = [30.0, 40.0, 55.0]; lit[:NADPH] = [50.0, 70.0, 65.0]
    metab_obs = Observable(fill(10.0, length(G._PANEL_METAB_ORDER)))
    fig = Figure(); ax = Axis(fig[1, 1])
    @test_nowarn G.draw_metab_panel!(ax, lit, metab_obs)
    @test ax.yscale[] === log10
    out = joinpath(tempdir(), "metab_panel_smoke.png"); @test_nowarn save(out, fig); @test isfile(out)
    @test_nowarn (metab_obs[] = fill(100.0, length(G._PANEL_METAB_ORDER)))
end
