using Documenter
using Literate
using GraphLab

generated_dir = joinpath(@__DIR__, "src", "generated")

literate_files = [
    "quickstart.jl",
    "bisection_methods.jl",
    "recursive_partitioning.jl",
]

for file in literate_files
    Literate.markdown(
        joinpath(@__DIR__, "literate", file),
        generated_dir;
        documenter=true,
    )
end

makedocs(
    sitename="GraphLab.jl",
    modules=[GraphLab],
    workdir=joinpath(@__DIR__, ".."),
    format=Documenter.HTML(
        repolink="https://github.com/lechekhabm/GraphLab.jl",
        collapselevel=1,
        prettyurls=true,
    ),
    repo="https://github.com/lechekhabm/GraphLab.jl",
    pages=[
        "Home" => "index.md",
        "Usage Guide" => "usage.md",
        "Examples" => [
            "Quick start" => "generated/quickstart.md",
            "Bisection methods" => "generated/bisection_methods.md",
            "Recursive partitioning" => "generated/recursive_partitioning.md",
        ],
        "API Reference" => "api.md",
        "Developers API Reference" => "dev_api.md",
    ],
    warnonly=true,
)

deploydocs(
    repo="github.com/lechekhabm/GraphLab.jl.git",
    branch="gh-pages",
    devbranch="main",
    versions=["stable", "v#.#.#", "dev"],
    forcepush=true,
    deploy_config=Documenter.GitHubActions(),
)