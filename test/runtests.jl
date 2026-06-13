using Test
@testset "GlycolysisPentosePathwayPhaseMap" begin
    include("test_grid_io.jl")
    include("test_flux_map.jl")
    include("test_metab_panel.jl")
    include("test_phasemap.jl")
    include("test_explorer.jl")
end
