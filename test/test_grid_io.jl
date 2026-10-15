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

@testset "nearest_cell (log-spaced)" begin
    xs = 10.0 .^ range(-3, 0, length = 4); ys = copy(xs)
    @test G.nearest_cell(xs, ys, 1e-3, 1e0) == (1, 4)
    @test G.nearest_cell(xs, ys, 10^(-1.6), 1e-3)[1] == 2
    @test G.nearest_cell(xs, ys, 10^(-1.4), 1e-3)[1] == 3
end

@testset "log_cell_edges" begin
    e = G.log_cell_edges(10.0 .^ range(-3, 0, length = 4))
    @test length(e) == 5
    @test issorted(e)
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
    phis = 10.0 .^ range(-3, 0, length = 12)
    df = DataFrame(i_nadph = repeat(1:12, outer = 12), i_r5p = repeat(1:12, inner = 12),
                   nadph_phi = repeat(phis, outer = 12), r5p_phi = repeat(phis, inner = 12))
    # nearest (log space) to NADPH φ 0.059, R5P φ 0.1068
    @test G.physiological_cell(df) == (7, 8)
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
