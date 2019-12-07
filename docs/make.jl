using Documenter, BeamFEA

makedocs(;
    modules=[BeamFEA],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/byuflowlab/BeamFEA.jl/blob/{commit}{path}#L{line}",
    sitename="BeamFEA.jl",
    authors="Andrew Ning <aning@byu.edu>",
    assets=String[],
)
