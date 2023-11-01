module FiveIGA

using Ferrite
using IGA
using Five

Ferrite.nnodes(::Type{<:IGA.BezierCell{shape,order,N}}) where {shape,order,N} = N

const IGALinearSolidElement{dim}   = LinearSolidElement{dim, <:BezierCellValues}
const IGAPart{dim,T}    = Part{dim,T,<:IGALinearSolidElement}

include("five_iga_geometry.jl")

Five._getquadraturerule(cv::BezierCellValues) = cv.cv_bezier.qr
Five._getinterpolation(cv::BezierCellValues) = cv.cv_bezier.ip

#default_geometry(::Part{2, Float64, LinearSolidElement{2, CellValues{VectorizedInterpolation{2, RefQuadrilateral, 2, IGAInterpolation{RefQuadrilateral, 2}}, Vec{2, Float64}, Tensor{2, 2, Float64, 4}, Tensor{2, 2, Float64, 4}, Float64, Vec{2, Float64}, QuadratureRule{RefQuadrilateral, Float64, 2}, IGAInterpolation{RefQuadrilateral, 2}}, PlaneStrain{2}}, LinearElastic}, ::BezierGrid{2, BezierCell{RefQuadrilateral, 2, 9}, Float64})

function Five.default_geometry(part::IGAPart, grid::BezierGrid)
    return IGASubGridGeometry(grid, part.cellset)
end

end # module  