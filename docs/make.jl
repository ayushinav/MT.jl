using MT
using Documenter

DocMeta.setdocmeta!(MT, :DocTestSetup, :(using MT); recursive=true)

makedocs(;
    modules=[MT],
    authors="Abhinav Pratap Singh",
    repo="https://github.com/ayushinav/MT.jl/blob/{commit}{path}#{line}",
    sitename="MT.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ayushinav.github.io/MT.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ayushinav/MT.jl",
    devbranch="main",
)
