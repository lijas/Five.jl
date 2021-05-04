using Documenter
using Five

makedocs(
    sitename = "Five",
    format = Documenter.HTML(),
    modules = [Five]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
