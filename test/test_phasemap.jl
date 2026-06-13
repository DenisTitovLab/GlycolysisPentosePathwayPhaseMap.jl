using Test, CairoMakie, GlycolysisPentosePathwayPhaseMap
CairoMakie.activate!()

@testset "phasemap() renders + saves" begin
    out = joinpath(tempdir(), "phasemap_smoke.png")
    fig = phasemap(save = out)            # uses the shipped demand_grid.csv
    @test fig isa Figure
    @test isfile(out)
end
