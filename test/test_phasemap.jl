using Test, CSV, DataFrames, CairoMakie, GlycolysisPentosePathwayPhaseMap
const G = GlycolysisPentosePathwayPhaseMap
CairoMakie.activate!()

@testset "phasemap() renders + saves" begin
    out = joinpath(tempdir(), "phasemap_smoke.png")
    fig = phasemap(save = out)            # uses the shipped demand_grid.csv
    @test fig isa Figure
    @test isfile(out)
    ax = only(filter(x -> x isa Axis && x.title[] == "Pentose-cycle index", fig.content))
    @test ax.xscale[] === identity && ax.yscale[] === identity     # equal-sized index cells
    @test ax.yticks[][2][end] == "100%"
end

@testset "phasemap default representatives" begin
    df = G.load_grid(G.default_demand_grid())
    reps = G._default_reps(df)
    @test length(reps) == 2
    @test occursin("physiological", reps[1][1])
    @test reps[1][2:3] == G.physiological_cell(df)
    rev = only(eachrow(df[(df.i_nadph .== reps[2][2]) .& (df.i_r5p .== reps[2][3]), :]))
    @test rev.retcode == :Terminated && rev.mode == "reverse"
end

@testset "phasemap() renders a grid without pool columns" begin
    df = CSV.read(G.default_demand_grid(), DataFrame)
    p = tempname() * ".csv"
    CSV.write(p, select(df, Not([:ADP_uM, :AMP_uM, :NAD_uM, :NADH_uM])))
    @test phasemap(grid = p) isa Figure
end
