using Test, CSV, CairoMakie, DataFrames, GlycolysisPentosePathwayPhaseMap
const G = GlycolysisPentosePathwayPhaseMap
CairoMakie.activate!()

@testset "load_metab_lit (shipped literature_pools.csv)" begin
    path = G.metab_lit_path()
    @test basename(path) == "literature_pools.csv"
    lit = G.load_metab_lit(path)
    @test Set(keys(lit)) == Set(G._PANEL_METAB_ORDER)
    @test all(v -> v isa Vector{Float64}, values(lit))
    @test all(v -> all(x -> isfinite(x) && x > 0, v), values(lit))
    @test all(s -> !isempty(lit[s]), filter(∉((:PGLn, :AMP)), collect(G._PANEL_METAB_ORDER)))
    raw = CSV.read(path, DataFrame)
    @test names(raw) == ["metabolite", "value_uM", "reference", "cell_type", "derived", "method"]
    @test Set(raw.method) ⊆ Set(["metabolomics", "sensor", "equilibrium"])
    left  = in(string.(G._FREE_EQ_TOTAL_PANEL)).(raw.metabolite)
    right = in(string.(G._FREE_POOL_PANEL)).(raw.metabolite)
    @test all(raw.method[left] .== "metabolomics") && !any(raw.derived[left])
    @test !any(raw.method[right] .== "metabolomics")
    @test all(raw.derived[raw.metabolite .== "NADP"]) && any(raw.metabolite .== "NADP")
    @test !any(raw.metabolite .== "AMP")
    @test all(s -> any(raw.metabolite .== string(s)), (:ADP, :NAD, :NADH, :NADPH))
    lo, hi = G._METAB_YLIMS                                   # no literature point is off-scale
    @test all(v -> lo <= v <= hi, raw.value_uM)
    @test isascii(read(path, String))
    @test !isfile(joinpath(G.DATA_DIR, "Metabolite_Concentrations.csv"))
    @test !isdefined(G, :_GAP_CORRECTION) && !isdefined(G, :_DHAP_CORRECTION)
end

@testset "load_metab_lit long format (synthetic)" begin
    path = tempname() * ".csv"
    CSV.write(path, DataFrame(metabolite = ["G6P", "G6P", "Foo"], value_uM = [1.0, 2.0, 3.0],
                              reference = ["a", "b", "c"], cell_type = ["x", "y", "z"],
                              derived = [false, false, false]))
    lit = G.load_metab_lit(path)
    @test lit[:G6P] == [1.0, 2.0]          # no conversion applied by the loader
    @test isempty(lit[:PGLn])
    @test !haskey(lit, :Foo)
end

@testset "draw_metab_panel! offscreen" begin
    lit = Dict(s => Float64[] for s in G._PANEL_METAB_ORDER)
    lit[:G6P] = [30.0, 40.0, 55.0]; lit[:NADPH] = [50.0, 70.0, 65.0]
    metab_obs = Observable(fill(10.0, length(G._PANEL_METAB_ORDER)))
    fig = Figure(); ax = Axis(fig[1, 1])
    @test_nowarn G.draw_metab_panel!(ax, lit, metab_obs)
    @test ax.yscale[] === log10
    @test ax.ylabel[] == G._METAB_YLABEL
    @test occursin("cytosolic water", G._METAB_YLABEL) && occursin("free", G._METAB_YLABEL)
    out = joinpath(tempdir(), "metab_panel_smoke.png"); @test_nowarn save(out, fig); @test isfile(out)
    @test_nowarn (metab_obs[] = fill(100.0, length(G._PANEL_METAB_ORDER)))
end
