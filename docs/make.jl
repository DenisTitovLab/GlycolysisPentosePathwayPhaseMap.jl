using Documenter, GlycolysisPentosePathwayPhaseMap
makedocs(
    sitename = "GlycolysisPentosePathwayPhaseMap.jl",
    modules  = [GlycolysisPentosePathwayPhaseMap],
    pages = ["Home" => "index.md", "Tutorial" => "tutorial.md",
             "The science" => "science.md", "Recompute" => "recompute.md",
             "API reference" => "api.md"],
    checkdocs = :exports,
)
deploydocs(repo = "github.com/DenisTitovLab/GlycolysisPentosePathwayPhaseMap.jl")
