module FiveIGA

using Ferrite
using IGA
using Five

Ferrite.nnodes(::Type{IGA.BezierCell{shape,order,N}}) where {shape,order,N} = N

#const IGAElement{dim}   = LinearSolidElement{dim, <:BezierCellValues}
#const IGAPart{dim,T}    = Part{dim,T,<:IGAElement}

include("five_iga_geometry.jl")

function Five.default_geometry(part, grid::BezierGrid)
    asdf
    return IGASubGridGeometry(grid, part.cellset)
end

end # module  