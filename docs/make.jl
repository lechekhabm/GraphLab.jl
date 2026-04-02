using Documenter
using Literate
using GraphLab  # Make sure it uses the locally developed version

# println("Defined functions: ", names(GraphPartitioning, all=true))
Literate.markdown(joinpath(@__DIR__, "literate", "ex1.jl"),
    joinpath(@__DIR__, "src", "generated");
    documenter=true)

Literate.markdown(joinpath(@__DIR__, "literate", "ex2.jl"),
    joinpath(@__DIR__, "src", "generated");
    documenter=true)

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
            "Example 1" => "generated/ex1.md",
            "Example 2" => "generated/ex2.md",
        ],
        "API Reference" => "api.md",
        "Developers API Reference" => "dev_api.md",
    ],
    warnonly=true,  # Prevent build failure due to missing docs
)


deploydocs(;
    repo="github.com/lechekhabm/GraphLab.jl.git",
    branch="gh-pages",
    devbranch="main",
    versions=["stable", "v#.#.#", "dev"],  # enable stable badge support
    forcepush=true,  # Ensure it force-pushes
    deploy_config=Documenter.GitHubActions()  # Adjust if your default branch is different
)