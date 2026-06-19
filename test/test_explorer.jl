using Test, CairoMakie, DataFrames, GlycolysisPentosePathwayPhaseMap
const G = GlycolysisPentosePathwayPhaseMap
CairoMakie.activate!()

@testset "log-space click → cell" begin
    xs = 10.0 .^ range(-3, 0, length = 4); ys = copy(xs)
    @test G.logclick_to_cell(xs, ys, (log10(xs[2]), log10(ys[3]))) == (2, 3)
    @test G.logclick_to_cell(xs, ys, (-1.4, -3.0))[1] == 3
end

@testset "build_explorer offscreen" begin
    df = G.load_grid(G.default_atp_grid())
    fig, selected, atp_level = G.build_explorer(df)
    @test selected[] == (cld(maximum(df.i_nadph), 2), cld(maximum(df.i_r5p), 2))
    atp_levels = sort(unique(df.atpase_frac))
    @test length(atp_levels) == 11
    @test atp_level[] == atp_levels[argmin(abs.(log10.(atp_levels) .- log10(0.10)))]
    out = joinpath(tempdir(), "explorer_smoke.png"); @test_nowarn save(out, fig); @test isfile(out)
    conv = first(filter(r -> r.retcode == :Terminated, eachrow(df)))
    @test_nowarn (selected[] = (conv.i_nadph, conv.i_r5p))
    @test_nowarn (atp_level[] = atp_levels[1])
    @test_nowarn save(out, fig)
    ax_metab = only(filter(x -> x isa Axis && x.ylabel[] == "[metabolite], µM", fig.content))
    for itr in (:rectanglezoom, :dragpan, :scrollzoom)
        @test !haskey(ax_metab.interactions, itr)
    end
end

@testset "explorer lookup indices" begin
    df = G.load_grid(G.default_atp_grid())
    nx = maximum(df.i_r5p); ny = maximum(df.i_nadph)

    # row index returns the same row the old per-interaction filter did
    idx = G._build_row_index(df)
    r = first(filter(x -> x.retcode == :Terminated, eachrow(df)))
    got = idx[(r.atpase_frac, r.i_nadph, r.i_r5p)]
    @test got.cycle_index == r.cycle_index
    @test got.r5p_phi == r.r5p_phi && got.nadph_phi == r.nadph_phi

    # zdict matches the old per-level zmatrix construction for EVERY ATP level
    zd = G._build_zdict(df, nx, ny)
    @test Set(keys(zd)) == Set(unique(df.atpase_frac))
    for atp in unique(df.atpase_frac)
        expected = fill(NaN, nx, ny)
        for rr in eachrow(df)
            (rr.atpase_frac == atp && rr.retcode == :Terminated) || continue
            expected[rr.i_r5p, rr.i_nadph] = rr.cycle_index
        end
        Z = zd[atp]
        @test size(Z) == (nx, ny)
        @test all((isnan.(Z) .& isnan.(expected)) .| (Z .== expected))
    end
end
