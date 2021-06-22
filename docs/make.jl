using Documenter
using Five

makedocs(
    sitename = "Five",
    format = Documenter.HTML(),
    modules = [Five]
)

deploydocs(
    repo = "github.com/lijas/Five.jl.git",
    push_preview=true,
    devbranch = "master"
)