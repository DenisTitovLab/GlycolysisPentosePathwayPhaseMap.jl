using Test, CairoMakie, DataFrames, GlycolysisPentosePathwayPhaseMap
const G = GlycolysisPentosePathwayPhaseMap
CairoMakie.activate!()

@testset "click → cell (equal-sized index cells)" begin
    @test G.click_to_cell(15, 18, (2.0, 3.0)) == (2, 3)
    @test G.click_to_cell(15, 18, (2.4, 3.6)) == (2, 4)
    @test G.click_to_cell(15, 18, (0.6, 18.4)) == (1, 18)     # edge clicks stay on the grid
    @test G.click_to_cell(15, 18, (-3.0, 40.0)) == (1, 18)
end

@testset "build_explorer offscreen" begin
    df = G.load_grid(G.default_atp_grid())
    fig, selected, atp_level = G.build_explorer(df)
    @test selected[] == G.physiological_cell(df)
    atp_levels = sort(unique(df.atpase_frac))
    @test length(atp_levels) == 11
    @test atp_level[] == atp_levels[argmin(abs.(log10.(atp_levels) .-
                                                log10(G._PHYSIOLOGICAL_DEMAND.atpase_frac)))]
    ax_heat = only(filter(x -> x isa Axis && occursin("Pentose-cycle", x.title[]), fig.content))
    # equal-sized cells: index positions, one percent-labelled tick per demand level
    @test ax_heat.xscale[] === identity && ax_heat.yscale[] === identity
    xs = sort(unique(df.r5p_phi)); ys = sort(unique(df.nadph_phi))
    @test ax_heat.xticks[] == (eachindex(xs), G.pct_label.(xs))
    @test ax_heat.yticks[] == (eachindex(ys), G.pct_label.(ys))
    out = joinpath(tempdir(), "explorer_smoke.png"); @test_nowarn save(out, fig); @test isfile(out)
    conv = first(filter(r -> r.retcode == :Terminated, eachrow(df)))
    ax_net = only(filter(x -> x isa Axis && occursin("ATP", x.title[]), fig.content))
    net_plots_before = copy(ax_net.scene.plots)
    @test_nowarn (selected[] = (conv.i_nadph, conv.i_r5p))
    @test_nowarn (atp_level[] = atp_levels[1])
    @test length(ax_net.scene.plots) == length(net_plots_before)
    @test all(a === b for (a, b) in zip(ax_net.scene.plots, net_plots_before))
    @test_nowarn save(out, fig)
    metab_axes = filter(x -> x isa Axis && x.yscale[] === log10 &&
                             (x.ylabel[] == G._METAB_YLABEL || occursin("free pool", x.title[])),
                        fig.content)
    @test length(metab_axes) == 2
    for ax in metab_axes, itr in (:rectanglezoom, :dragpan, :scrollzoom)
        @test !haskey(ax.interactions, itr)
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
