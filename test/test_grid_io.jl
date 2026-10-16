using Test, CSV, DataFrames, GlycolysisPentosePathwayPhaseMap
const G = GlycolysisPentosePathwayPhaseMap

@testset "pentose_cycle_index / classify_mode" begin
    @test G.pentose_cycle_index((; V_G6PD = 1.0, V_R5Pase = 0.0)) ≈ 1.0
    @test G.pentose_cycle_index((; V_G6PD = 1.0, V_R5Pase = 1.5)) ≈ -0.5
    @test isnan(G.pentose_cycle_index((; V_G6PD = 0.0, V_R5Pase = 0.1)))
    @test G.classify_mode((; V_G6PD = 1.0, V_R5Pase = 0.1)) == :cycle
    @test G.classify_mode((; V_G6PD = 1.0, V_R5Pase = 0.7)) == :linear
    @test G.classify_mode((; V_G6PD = 1.0, V_R5Pase = 1.4)) == :reverse
    @test G.classify_mode((; V_G6PD = 0.0, V_R5Pase = 0.1)) == :undetermined
end

@testset "default demand axes" begin
    @test G._DEFAULT_R5P_PHIS == [0.001, 0.003, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.07, 0.10,
                                  0.15, 0.20, 0.30, 0.50, 1.0]
    @test G._DEFAULT_NADPH_PHIS == sort([G._DEFAULT_R5P_PHIS; 0.6; 0.7; 0.8])
    @test length(G._DEFAULT_NADPH_PHIS) == 18
    for phis in (G._DEFAULT_R5P_PHIS, G._DEFAULT_NADPH_PHIS)
        @test issorted(phis) && allunique(phis)
        @test count(<(0.01), phis) == 2                  # sparse below 1 %
        @test count(x -> 0.01 <= x <= 0.2, phis) == 10   # dense over 1-20 %
    end
end

@testset "nearest_cell (log-spaced)" begin
    xs = 10.0 .^ range(-3, 0, length = 4); ys = copy(xs)
    @test G.nearest_cell(xs, ys, 1e-3, 1e0) == (1, 4)
    @test G.nearest_cell(xs, ys, 10^(-1.6), 1e-3)[1] == 2
    @test G.nearest_cell(xs, ys, 10^(-1.4), 1e-3)[1] == 3
    ys = G._DEFAULT_NADPH_PHIS                       # snaps onto the added high-NADPH rows
    @test G.nearest_cell(xs, ys, 1e-3, 0.68)[2] == 16
end

@testset "index_position" begin
    lv = [0.01, 0.1, 1.0]
    @test G.index_position(lv, 0.1) == 2.0
    @test G.index_position(lv, sqrt(0.1)) ≈ 2.5            # log interpolation, not linear
    @test G.index_position(lv, 1e-4) == 1.0 && G.index_position(lv, 5.0) == 3.0   # clamped
    # physiological demand lands just past its nearest cell (0.05, 0.10) on the default axes
    @test 8 < G.index_position(G._DEFAULT_NADPH_PHIS, G._PHYSIOLOGICAL_DEMAND.nadph_phi) < 8.5
    @test 10 < G.index_position(G._DEFAULT_R5P_PHIS, G._PHYSIOLOGICAL_DEMAND.r5p_phi) < 10.5
end

@testset "pct_label" begin
    @test G.pct_label.(G._DEFAULT_R5P_PHIS) ==
          ["0.1%", "0.3%", "1%", "1.5%", "2%", "3%", "4%", "5%", "7%", "10%", "15%", "20%",
           "30%", "50%", "100%"]
    @test G.pct_label(0.0123) == "1.23%"
end

@testset "row_to_pools ordering + NaN passthrough" begin
    order = G._PANEL_METAB_ORDER
    cols = Dict(Symbol(string(s), "_uM") => Float64(i) for (i, s) in enumerate(order))
    row = first(eachrow(DataFrame(cols)))
    @test G.row_to_pools(row) == Float64.(1:length(order))
    cols_nan = Dict(Symbol(string(s), "_uM") => NaN for s in order)
    @test all(isnan, G.row_to_pools(first(eachrow(DataFrame(cols_nan)))))
end

@testset "reversed R5P sink -> undetermined" begin
    # index > 1 needs V_R5Pase < 0: the sink is supplying R5P (ATP collapse), not recycling
    @test isnan(G.pentose_cycle_index((; V_G6PD = 1.0, V_R5Pase = -0.5)))
    @test G.classify_mode((; V_G6PD = 1.0, V_R5Pase = -0.5)) === :undetermined
    @test G.pentose_cycle_index((; V_G6PD = 1.0, V_R5Pase = 0.0)) ≈ 1.0
end

@testset "physiological_cell" begin
    xs = G._DEFAULT_R5P_PHIS; ys = G._DEFAULT_NADPH_PHIS
    df = DataFrame(i_nadph = repeat(eachindex(ys), outer = length(xs)),
                   i_r5p = repeat(eachindex(xs), inner = length(ys)),
                   nadph_phi = repeat(ys, outer = length(xs)),
                   r5p_phi = repeat(xs, inner = length(ys)))
    # nearest (log space) to NADPH φ 0.059, R5P φ 0.1068: the (0.05, 0.10) cell
    @test G.physiological_cell(df) == (8, 10)
    @test (ys[8], xs[10]) == (0.05, 0.10)
end

@testset "shipped grids carry the viewer's columns" begin
    for path in (G.default_atp_grid(), G.default_demand_grid())
        cols = names(CSV.read(path, DataFrame; limit = 1))
        for k in G._PPP_FLUX_KEYS
            @test "$(k)_uM_min" in cols
        end
        for s in G._PANEL_METAB_ORDER
            @test "$(s)_uM" in cols
        end
        @test "retcode" in cols && "cycle_index" in cols && "mode" in cols
    end
end

@testset "shipped grids carry the cofactor pools and the F26BP setting" begin
    for path in (G.default_atp_grid(), G.default_demand_grid())
        cols = names(CSV.read(path, DataFrame; limit = 1))
        for s in (:ADP, :AMP, :NAD, :NADH)
            @test "$(s)_uM" in cols
        end
    end
    @test G.read_provenance()["grids"]["f26bp_uM"] == G._PHYSIOLOGICAL_DEMAND.f26bp_uM
end

@testset "shipped provenance matches shipped data" begin
    p = G.read_provenance()
    @test p["grids"]["ppp_commit"] == "fcccae8a196ce9cfba23fae91ca77c4be4540807"
    @test p["grids"]["ppp_dirty"] == false
    @test p["grids"]["atp_grid_rows"] == nrow(G.load_grid(G.default_atp_grid()))
    @test p["grids"]["demand_grid_rows"] == nrow(G.load_grid(G.default_demand_grid()))
    @test p["grids"]["atp_grid_nonconverged"] ==
          count(!=(:Terminated), G.load_grid(G.default_atp_grid()).retcode)
    @test p["grids"]["atpase_frac"] == G._PHYSIOLOGICAL_DEMAND.atpase_frac
    @test p["grids"]["atp_fracs"][end] ≈ G._PHYSIOLOGICAL_DEMAND.atpase_frac
    @test !isempty(p["literature"]["ppp_experiments_commit"])
    # the shipped grids sit on the default linear demand axes
    @test p["grids"]["nadph_phis"] == G._DEFAULT_NADPH_PHIS
    @test p["grids"]["r5p_phis"] == G._DEFAULT_R5P_PHIS
    for path in (G.default_atp_grid(), G.default_demand_grid())
        df = G.load_grid(path)
        @test sort(unique(df.nadph_phi)) == G._DEFAULT_NADPH_PHIS
        @test sort(unique(df.r5p_phi)) == G._DEFAULT_R5P_PHIS
    end
end

@testset "shipped grids: sinks within capacity, no index > 1, physiological cell alive" begin
    for path in (G.default_atp_grid(), G.default_demand_grid())
        df = G.load_grid(path)
        ok = df.retcode .== :Terminated
        # fluxes cached in µM/min and Vmax columns in M/min: a sink never exceeds its capacity
        @test all(df.V_R5Pase_uM_min[ok] .<= df.R5Pase_Vmax[ok] .* 1e6 .* (1 + 1e-6))
        @test all(df.V_NADPHox_uM_min[ok] .<= df.NADPHox_Vmax[ok] .* 1e6 .* (1 + 1e-6))
        @test !any(x -> x > 1 + 1e-9, filter(!isnan, df.cycle_index))
    end
    d2 = G.load_grid(G.default_demand_grid())
    in_, ir = G.physiological_cell(d2)
    row = only(eachrow(d2[(d2.i_nadph .== in_) .& (d2.i_r5p .== ir), :]))
    @test row.retcode == :Terminated && row.ATP_uM > 1000
end

@testset "panel tuples partition the panel order" begin
    l, r = G._FREE_EQ_TOTAL_PANEL, G._FREE_POOL_PANEL
    @test r == (:ADP, :AMP, :NAD, :NADH, :NADP, :NADPH)
    @test length(l) == 14 && isempty(intersect(l, r))
    @test G._PANEL_METAB_ORDER == (l..., r...)
    for path in (G.default_atp_grid(), G.default_demand_grid())
        cols = names(CSV.read(path, DataFrame; limit = 1))
        @test all("$(s)_uM" in cols for s in G._PANEL_METAB_ORDER)
    end
end

@testset "grid column checks" begin
    df = CSV.read(G.default_atp_grid(), DataFrame)
    p = tempname() * ".csv"
    CSV.write(p, select(df, Not(:V_G6PD_uM_min)))
    err = try G.load_grid(p); nothing catch e; e end
    @test err isa ArgumentError && occursin("V_G6PD_uM_min", err.msg)

    p2 = tempname() * ".csv"
    CSV.write(p2, select(df, Not([:ADP_uM, :NADH_uM])))
    old = G.load_grid(p2)                                   # load_grid does not need pools
    err2 = try G.build_explorer(old); nothing catch e; e end
    @test err2 isa ArgumentError
    @test occursin("ADP_uM", err2.msg) && occursin("NADH_uM", err2.msg)
    @test occursin("regenerate", err2.msg)

    p3 = tempname() * ".csv"
    CSV.write(p3, select(df, Not(:atpase_frac)))
    err3 = try G.build_explorer(G.load_grid(p3)); nothing catch e; e end
    @test err3 isa ArgumentError && occursin("atpase_frac", err3.msg)
end
