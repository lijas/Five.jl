using Documenter
using Five

include("generate.jl")

GENERATEDEXAMPLES = [joinpath("examples", f) for f in (
    "bar_example.md",
    "beam_example.md",
    "enf_example.md"
    )]

makedocs(
    sitename = "Five",
    format = Documenter.HTML(),
    #modules = [Five],
    pages = Any[
        "Home" => "index.md",
        "essentials.md",
        "parts.md",
        "Elements"  => [
            "elements/elements_overview.md",
            "elements/solid_element.md"
            ],
        "Solvers" => [
            "solvers/solver_overview.md",
            "solvers/crisfield_solver.md",
            "solvers/local_dissipation_solver.md",
            ],
        "Examples" => GENERATEDEXAMPLES,
    ]
)

deploydocs(
    repo = "github.com/lijas/Five.jl.git",
    push_preview=true,
    devbranch = "master"
)