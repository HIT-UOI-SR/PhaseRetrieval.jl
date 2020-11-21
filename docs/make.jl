using PhaseRetrieval
using Documenter

makedocs(;
    modules=[PhaseRetrieval],
    authors="Yong-an Lu <miroox@outlook.com>",
    repo="https://github.com/HIT-UOI-SR/PhaseRetrieval.jl/blob/{commit}{path}#L{line}",
    sitename="PhaseRetrieval.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://HIT-UOI-SR.github.io/PhaseRetrieval.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/HIT-UOI-SR/PhaseRetrieval.jl",
)
