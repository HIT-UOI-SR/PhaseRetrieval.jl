using PhaseRetrieval
using Documenter

makedocs(;
    modules=[PhaseRetrieval],
    authors="Yong-an Lu <miroox@outlook.com>",
    repo=Documenter.Remotes.GitHub("HIT-UOI-SR", "PhaseRetrieval.jl"),
    sitename="PhaseRetrieval.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://HIT-UOI-SR.github.io/PhaseRetrieval.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
