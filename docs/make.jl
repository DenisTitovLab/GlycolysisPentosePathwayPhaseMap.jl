using GlycolysisPentosePathwayPhaseMap
using Documenter

DocMeta.setdocmeta!(GlycolysisPentosePathwayPhaseMap, :DocTestSetup, :(using GlycolysisPentosePathwayPhaseMap); recursive=true)

makedocs(;
    modules=[GlycolysisPentosePathwayPhaseMap],
    authors="James Mbata",
    sitename="GlycolysisPentosePathwayPhaseMap.jl",
    format=Documenter.HTML(;
        canonical="https://DenisTitovLab.github.io/GlycolysisPentosePathwayPhaseMap.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/DenisTitovLab/GlycolysisPentosePathwayPhaseMap.jl",
    devbranch="main",
)
