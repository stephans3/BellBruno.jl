using BellBruno
using Documenter

DocMeta.setdocmeta!(BellBruno, :DocTestSetup, :(using BellBruno); recursive=true)

makedocs(;
    modules=[BellBruno],
    authors="Stephan Scholz",
    repo="https://github.com/stephans3/BellBruno.jl/blob/{commit}{path}#{line}",
    sitename="BellBruno.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://stephans3.github.io/BellBruno.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/stephans3/BellBruno.jl",
    devbranch="main",
)
