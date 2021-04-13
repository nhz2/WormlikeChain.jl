using WormlikeChain
using Documenter

DocMeta.setdocmeta!(WormlikeChain, :DocTestSetup, :(using WormlikeChain); recursive=true)

makedocs(;
    modules=[WormlikeChain],
    authors="Nathan Zimmerberg <nhz2@cornell.edu> and contributors",
    repo="https://github.com/nhz2/WormlikeChain.jl/blob/{commit}{path}#{line}",
    sitename="WormlikeChain.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://nhz2.github.io/WormlikeChain.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/nhz2/WormlikeChain.jl",
    devbranch = "main"
)
