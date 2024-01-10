module FiveIGA

using Ferrite
using IGA
using Five

Ferrite.nnodes(::Type{<:IGA.BezierCell{Ferrite.RefHypercube{dim},order}}) where {dim,order} = (order+1)^dim

const IGALinearSolidElement{dim}   = LinearSolidElement{dim, <:BezierCellValues}
const IGASolidElement{dim}         = SolidElement{dim, <:BezierCellValues}
const IGAElement{dim}             = Union{IGALinearSolidElement{dim}, IGASolidElement{dim}}

const IGAPart{dim,T}    = Part{dim,T,<:IGAElement}


include("five_iga_geometry.jl")

Five._getquadraturerule(cv::BezierCellValues) = cv.cv_bezier.qr
Five._getinterpolation(cv::BezierCellValues) = cv.cv_bezier.ip

#default_geometry(::Part{2, Float64, LinearSolidElement{2, CellValues{VectorizedInterpolation{2, RefQuadrilateral, 2, IGAInterpolation{RefQuadrilateral, 2}}, Vec{2, Float64}, Tensor{2, 2, Float64, 4}, Tensor{2, 2, Float64, 4}, Float64, Vec{2, Float64}, QuadratureRule{RefQuadrilateral, Float64, 2}, IGAInterpolation{RefQuadrilateral, 2}}, PlaneStrain{2}}, LinearElastic}, ::BezierGrid{2, BezierCell{RefQuadrilateral, 2, 9}, Float64})

function IGAPart(; 
    material::M,
    element::E,
    cellset,
    ) where {E<:IGAElement,M}

    dim = Ferrite.getdim(element)
    T = Float64

    _set = collect(cellset)
    sort!(_set) # YOLO

    coords_t = IGA.BezierCoords{dim,T} 
    cache_t = Five.PartCache{dim,T,E,coords_t}

    return Part{dim,T,E,M,cache_t}(
        material, 
        _set,
        Vector{Int}[],
        element, 
        cache_t[])
end

function Five.create_part_cache(element::E) where E<:IGAElement
    dim = Ferrite.getdim(element)
    T = Float64
    _ndofs   = ndofs(element)
    _nnodes  = Ferrite.nnodes(getcelltype(element))

    coords_t = IGA.BezierCoords{dim,T} 
    return Five.PartCache{dim,T,E,coords_t}(
        zeros(T,_ndofs), 
        zeros(T,_ndofs), 
        zeros(T,_ndofs), 
        zeros(T,_ndofs), 
        zeros(T,_ndofs,_ndofs), 
        zeros(Int,_ndofs), 
        IGA.zero_bezier_coord(dim,T,_nnodes),
        deepcopy(element)
    )
end

function Five.default_geometry(part::IGAPart, grid::BezierGrid)
    return IGASubGridGeometry(grid, part.cellset)
end

end # module  