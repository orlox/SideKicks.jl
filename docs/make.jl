using SideKicks
using Documenter
using Literate

# Parse examples using Literate
pkg_path = pkgdir(SideKicks)
@show

function ignore_code_blocks(content)
    content = replace(content, "##\n" => "\n")  # remove code blocks
    content = replace(content, "###" => "##")  # make level 3 headers level 2
end

Literate.markdown(pkg_path * "/examples/run_inference_vfts243.jl", pkg_path * "/docs/src/", preprocess=ignore_code_blocks, name="1_run_inference")
Literate.markdown(pkg_path * "/examples/plot_results_vfts243.jl", pkg_path * "/docs/src/", preprocess=ignore_code_blocks, name="2_plot_results")

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
        "Examples" => ["1_run_inference.md", "2_plot_results.md"]
    ],
)

deploydocs(;
    repo="github.com/orlox/SideKicks.jl",
    devbranch="main",
)
