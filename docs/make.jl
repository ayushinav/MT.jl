using Pkg
Pkg.activate("../MT/.")

using MT
using Documenter

DocMeta.setdocmeta!(MT, :DocTestSetup, :(using MT); recursive=true)

makedocs(;
    modules=[MT],
    authors="Abhinav Pratap Singh",
    repo="https://github.com/ayushinav/MT", #.jl/blob/{commit}{path}#{line}",
    sitename="MT.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ayushinav.github.io/MT",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
        "forward" => "forward.md",
        "inverse" => "inverse.md",
        "model" => "model.md",
    ],
)

deploydocs(;
    repo="github.com/ayushinav/MT.git",
    devbranch="main",
)
