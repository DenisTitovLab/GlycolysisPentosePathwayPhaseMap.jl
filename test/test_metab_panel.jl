using Test, CSV, CairoMakie, DataFrames, GlycolysisPentosePathwayPhaseMap
const G = GlycolysisPentosePathwayPhaseMap
CairoMakie.activate!()

@testset "load_metab_lit (shipped literature_pools.csv)" begin
    path = G.metab_lit_path()
    @test basename(path) == "literature_pools.csv"
    lit = G.load_metab_lit(path)
    @test Set(keys(lit)) == Set(G._PANEL_METAB_ORDER)
    @test all(v -> v isa Vector{G._LitPoint}, values(lit))
    @test all(v -> all(p -> isfinite(p.value) && p.value > 0, v), values(lit))
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
    @test all(p -> p.method === :metabolomics, lit[:G6P])
    @test only(lit[:NADP]).derived
    @test all(s -> any(raw.metabolite .== string(s)), (:ADP, :NAD, :NADH, :NADPH))
    lo, hi = G._METAB_YLIMS                                   # no literature point is off-scale
    @test all(v -> lo <= v <= hi, raw.value_uM)
    @test isascii(read(path, String))
    @test !isfile(joinpath(G.DATA_DIR, "Metabolite_Concentrations.csv"))
    @test !isdefined(G, :_GAP_CORRECTION) && !isdefined(G, :_DHAP_CORRECTION)
end

@testset "load_metab_lit long format (synthetic)" begin
    path = tempname() * ".csv"
    CSV.write(path, DataFrame(metabolite = ["G6P", "G6P", "NADH", "Foo"], value_uM = [1.0, 2.0, 0.3, 3.0],
                              reference = ["a", "b", "c", "d"], cell_type = ["x", "y", "z", "w"],
                              derived = [false, false, true, false],
                              method = ["metabolomics", "metabolomics", "equilibrium", "sensor"]))
    lit = G.load_metab_lit(path)
    @test [p.value for p in lit[:G6P]] == [1.0, 2.0]   # no conversion applied by the loader
    @test only(lit[:NADH]) == (value = 0.3, derived = true, method = :equilibrium)
    @test isempty(lit[:PGLn]) && isempty(lit[:AMP])
    @test !haskey(lit, :Foo)

    wrong(m, meth) = (p = tempname() * ".csv";
        CSV.write(p, DataFrame(metabolite = [m], value_uM = [1.0], reference = ["a"],
                               cell_type = ["x"], derived = [false], method = [meth])); p)
    @test_throws ArgumentError G.load_metab_lit(wrong("G6P", "sensor"))       # left, not metabolomics
    @test_throws ArgumentError G.load_metab_lit(wrong("NADH", "metabolomics")) # total pool on the right
    @test_throws ArgumentError G.load_metab_lit(wrong("G6P", "guess"))        # unknown method
end

@testset "_free_pool_point_style (basic-panel marker key)" begin
    s = G._free_pool_point_style(false, :sensor)
    @test s.marker == :circle && s.color == :black
    @test G._free_pool_point_style(true, :sensor).color == :transparent
    @test G._free_pool_point_style(false, :equilibrium).marker == :diamond
end

@testset "_offscale_marks" begin
    m = G._offscale_marks([10.0, 1e-6, 0.0, -1e-18, NaN, 1e7]; lims = (1e-4, 1e5))
    @test m.on == [Point2f(1, 10.0)]
    @test [p[1] for p in m.below] == [2, 3, 4]               # tiny, zero and negative: floor
    @test [p[1] for p in m.above] == [6]                     # NaN draws nothing
    @test length(m.labels) == length(m.label_pos) == 4
    @test m.labels[2] == "0.0e+00"
    @test m.labels[1] == "1.0e-06"                           # 1e-6 entry, not rounded onto the floor
    @test all(p -> p[2] ≈ 1e-4 * 1.5, m.below) && all(p -> p[2] ≈ 1e5 / 1.5, m.above)
end

@testset "draw_metab_split! offscreen" begin
    lit = Dict(s => G._LitPoint[] for s in G._PANEL_METAB_ORDER)
    lit[:G6P] = [(value = v, derived = false, method = :metabolomics) for v in (30.0, 40.0, 55.0)]
    lit[:NADH] = [(value = 0.13, derived = false, method = :sensor),
                  (value = 0.34, derived = true, method = :equilibrium)]
    n = length(G._PANEL_METAB_ORDER)
    pools = fill(10.0, n); pools[2] = 1e-7; pools[16] = 0.0; pools[17] = NaN   # G6P off, AMP 0, NAD NaN
    metab_obs = Observable(pools)
    fig = Figure(size = (1500, 400)); gl = GridLayout(fig[1, 1])
    axes = G.draw_metab_split!(gl, lit, metab_obs)
    @test axes.ax_left.yscale[] === log10 && axes.ax_right.yscale[] === log10
    @test axes.ax_left.ylabel[] == G._METAB_YLABEL == "[metabolite], µM (cytosolic water)"
    @test occursin("no data", axes.ax_right.xticks[][2][2])                   # AMP column
    out = joinpath(tempdir(), "metab_split_smoke.png"); @test_nowarn save(out, fig); @test isfile(out)

    nplots = (length(axes.ax_left.scene.plots), length(axes.ax_right.scene.plots))
    for k in (3, 0, 5)                                        # off-scale count changes per update
        p = fill(10.0, n); p[1:k] .= 1e-9; p[end-k+1:end] .= 1e9
        @test_nowarn (metab_obs[] = p)
        @test_nowarn save(out, fig)
        @test (length(axes.ax_left.scene.plots), length(axes.ax_right.scene.plots)) == nplots
    end
    @test axes.ax_right.finallimits[].origin[1] ≈ 0.5
    @test axes.ax_right.finallimits[].widths[1] ≈ length(G._FREE_POOL_PANEL) + 0.3
    @test axes.ax_left.finallimits[].widths[1] ≈ length(G._FREE_EQ_TOTAL_PANEL) + 0.3
end
