using Documenter
using Five

makedocs(
    sitename = "Five",
    format = Documenter.HTML(),
    #modules = [Five],
    pages = Any[
        "Home" => "index.md",
        "essentials.md",
        "parts.md",
        "elements.md",
        "Solvers" => [
            "solvers/crisfield_solver.md",
            "solvers/local_dissipation_solver.md",
            ],
    ]
)

deploydocs(
    repo = "github.com/lijas/Five.jl.git",
    push_preview=true,
    devbranch = "master"
)