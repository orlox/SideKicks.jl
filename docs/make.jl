using SideKicks
using Documenter

DocMeta.setdocmeta!(SideKicks, :DocTestSetup, :(using SideKicks); recursive=true)

makedocs(;
    modules=[SideKicks],
    authors="Pablo Marchant <pamarca@gmail.com> and contributors",
    sitename="SideKicks.jl",
    format=Documenter.HTML(;
        canonical="https://orlox.github.io/SideKicks.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/orlox/SideKicks.jl",
    devbranch="main",
)
