export CohesiveCell

"""
"""

_get_interface_cell_shape(::Type{RefLine}) = RefQuadrilateral
_get_interface_cell_shape(::Type{RefTriangle}) = RefPrism
_get_interface_cell_shape(::Type{RefQuadrilateral}) = RefHexahedron

"""
CohesiveZoneInterpolation{RefLine,2}
3-----6-----4
|           |
1-----5-----2

refshape:       RefLine            #Not sure if it should be considered a RefLine or RefQuadrilateral
order:          2
ref_cell_shape: RefQuadrilateral   #This is needed for Ferrite.InterpolationInfo

"""
struct CohesiveZoneInterpolation{refshape,order,I<:Ferrite.ScalarInterpolation} <: Ferrite.ScalarInterpolation{refshape,order} 
    interpolation::I
    function CohesiveZoneInterpolation(interpolation)
        refshape_ip = Ferrite.getrefshape(interpolation)
        ref_cell_shape = _get_interface_cell_shape(refshape_ip)
        order = Ferrite.getorder(interpolation)
        new{ref_cell_shape,order,typeof(interpolation)}(interpolation)
    end
end

Ferrite.nfaces(::CohesiveZoneInterpolation{RefQuadrilateral}) = 2
Ferrite.nedges(::CohesiveZoneInterpolation{RefQuadrilateral}) = 0
Ferrite.nvertices(::CohesiveZoneInterpolation{RefQuadrilateral}) = 4

Ferrite.getnbasefunctions(ip::CohesiveZoneInterpolation) = getnbasefunctions(ip.interpolation)*2
Ferrite.vertexdof_indices(::CohesiveZoneInterpolation{RefQuadrilateral}) = ((1,),(2,),(3,),(4,))
Ferrite.facedof_indices(::CohesiveZoneInterpolation{RefQuadrilateral, 1}) = ((1,2), (3,4))
Ferrite.facedof_indices(::CohesiveZoneInterpolation{RefQuadrilateral, 2}) = ((1,2,5), (3,4,6))
Ferrite.facedof_interior_indices(::CohesiveZoneInterpolation{RefQuadrilateral, 2}) = ((5,),(6,))

Ferrite.adjust_dofs_during_distribution(::CohesiveZoneInterpolation{RefQuadrilateral,1}) = false
Ferrite.adjust_dofs_during_distribution(::CohesiveZoneInterpolation{RefQuadrilateral,2}) = false

get_cz_refshape(ip::CohesiveZoneInterpolation) = Ferrite.getrefshape(ip.interpolation)

_mapper(::Lagrange{RefLine,1}, _i::Int) = [1,2,1,2][_i]
_mapper(::Lagrange{RefLine,2}, _i::Int) = [1,2,1,2,3,3][_i]
_mapper(::Lagrange{RefQuadrilateral,1}, _i::Int) = [1,2,3,4,1,2,3,4][_i]
_mapper(::Lagrange{RefTriangle,1}, _i::Int) = [1,2,3,1,2,3][_i]
function mid_surf_value(ip::CohesiveZoneInterpolation, ξ::Vec{dim_p}, _i::Int) where dim_p
    _i<=getnbasefunctions(ip) || throw(ArgumentError("no shape function $_i for interpolation $ip"))
    i = _mapper(ip.interpolation, _i)

    return Ferrite.shape_value(ip.interpolation, ξ, i)*0.5
end

_sign_mapper(::Lagrange{RefLine,1}, _i::Int) = [-1,-1,+1,+1][_i]
_sign_mapper(::Lagrange{RefLine,2}, _i::Int) = [-1,-1,+1,+1,-1,+1][_i]
_sign_mapper(::Lagrange{RefQuadrilateral,1}, _i::Int) = [-1,-1,-1,-1,+1,+1,+1,+1][_i]
_sign_mapper(::Lagrange{RefTriangle,1}, _i::Int) = [-1,-1,-1,+1,+1,+1][_i]
function Ferrite.shape_value(ip::CohesiveZoneInterpolation, ξ::Vec{dim_p}, _i::Int) where dim_p
    _i<=getnbasefunctions(ip) || throw(ArgumentError("no shape function $_i for interpolation $ip"))
    sign = _sign_mapper(ip.interpolation, _i)
    i = _mapper(ip.interpolation, _i)
    return sign*Ferrite.shape_value(ip.interpolation, ξ, i)
end

Ferrite._mass_qr(::CohesiveZoneInterpolation{RefLine,1}) = QuadratureRule{RefLine}(1)

function Ferrite.reference_coordinates(::CohesiveZoneInterpolation{RefLine,1})
    return [Vec{1,Float64}((-1.0,)),
            Vec{1,Float64}((1.0,))]
end

"""

"""

abstract type AbstractCohesiveCell{refshape} <: Ferrite.AbstractCell{refshape} end
struct CZQuadrilateral          <: AbstractCohesiveCell{RefQuadrilateral} nodes::NTuple{4,Int} end
struct CZQuadraticQuadrilateral <: AbstractCohesiveCell{RefQuadrilateral} nodes::NTuple{6,Int} end
#struct CZTriangle      <: AbstractCohesiveCell{RefTriangle} nodes::NTuple{6,Int} end

Ferrite.vertices(c::AbstractCohesiveCell{RefQuadrilateral})= (c.nodes[1], c.nodes[2], c.nodes[3], c.nodes[4])
Ferrite.faces(c::AbstractCohesiveCell{RefQuadrilateral})   = ((c.nodes[1], c.nodes[2]), (c.nodes[3], c.nodes[4])) 

#Ferrite.vertices(c::CZQuadrilateral) = (c.nodes[1], c.nodes[2], c.nodes[3], c.nodes[4], c.nodes[5], c.nodes[6], c.nodes[7], c.nodes[8])
#Ferrite.faces(c::CZQuadrilateral) = ((c.nodes[1], c.nodes[2], c.nodes[3], c.nodes[4]), (c.nodes[5], c.nodes[6], c.nodes[7], c.nodes[8])) 
#Ferrite.vertices(c::CZTriangle) = (c.nodes[1], c.nodes[2], c.nodes[3], c.nodes[4], c.nodes[5], c.nodes[6],)
#Ferrite.faces(c::CZTriangle) = ((c.nodes[1], c.nodes[2], c.nodes[3], ), (c.nodes[4], c.nodes[6], c.nodes[6],)) 


Ferrite.default_interpolation(::Type{CZQuadrilateral}) = CohesiveZoneInterpolation(Lagrange{RefLine,1}())
Ferrite.default_interpolation(::Type{CZQuadraticQuadrilateral}) = CohesiveZoneInterpolation(Lagrange{RefLine,2}())
#Ferrite.default_interpolation(::Type{CZQuadrilateral}) = CohesiveZoneInterpolation(Lagrange{RefQuadrilateral,1}())
#Ferrite.default_interpolation(::Type{CZTriangle}) = CohesiveZoneInterpolation(Lagrange{RefTriangle,1}())

Ferrite.cell_to_vtkcell(::Type{CZQuadrilateral}) = Ferrite.VTKCellTypes.VTK_QUAD
Ferrite.cell_to_vtkcell(::Type{CZQuadraticQuadrilateral}) = Ferrite.VTKCellTypes.VTK_QUAD
#errite.cell_to_vtkcell(::Type{CZQuadrilateral}) = Ferrite.VTKCellTypes.VTK_HEXAHEDRON
#Ferrite.cell_to_vtkcell(::Type{CZTriangle}) = Ferrite.VTKCellTypes.VTK_WEDGE

Ferrite.nodes_to_vtkorder(cell::CZQuadraticQuadrilateral) = cell.nodes[[1,2,4,3]] #dont use center nodes?

"""
    SurfaceVectorValues

Contains shape values for cohesive zone elements
"""
struct SurfaceVectorValues{dim_p,dim_s,T<:Real,M2,refshape<:Ferrite.AbstractRefShape,ip} <: Ferrite.AbstractCellValues
    N::Matrix{Vec{dim_s,T}}
    dNdξ::Matrix{Tensor{2,dim_s,T,M2}}
    dNdx::Matrix{Tensor{2,dim_s,T,M2}}
    R::Vector{Tensor{2,dim_s,T,M2}}
    detJdA::Vector{T}
    M::Matrix{T}  # Shape values for geometric interp
    dMdξ::Matrix{Vec{dim_p,T}}
    qr::QuadratureRule{refshape,T}
    covar_base::Vector{Tensor{2,dim_s,T}}
    ip::Interpolation
end

Ferrite.getnquadpoints(cv::SurfaceVectorValues) = length(cv.qr.weights)

function SurfaceVectorValues(quad_rule::QuadratureRule, func_interpol::CohesiveZoneInterpolation, geom_interpol::CohesiveZoneInterpolation=func_interpol)
    SurfaceVectorValues(Float64, quad_rule, func_interpol, geom_interpol)
end

@inline getR(cv::SurfaceVectorValues, qp::Int) = cv.R[qp]

function SurfaceVectorValues(::Type{T},
    quad_rule::QuadratureRule{refshape},
    func_interpol::CohesiveZoneInterpolation,
    geom_interpol::CohesiveZoneInterpolation=func_interpol) where {T,refshape}

    @assert Ferrite.getdim(func_interpol) == Ferrite.getdim(geom_interpol)
    @assert Ferrite.getrefshape(func_interpol) == Ferrite.getrefshape(geom_interpol)
    dim_p = Ferrite.getdim(func_interpol)-1
    dim_s = dim_p+1

    n_qpoints = length(Ferrite.getweights(quad_rule))

    # Function interpolation
    n_func_basefuncs = getnbasefunctions(func_interpol) * dim_s #Note, multipy with two
    N    = fill(zero(Vec{dim_s,T}) * T(NaN), n_func_basefuncs, n_qpoints)
    dNdξ = fill(zero(Tensor{2,dim_s,T}) * T(NaN), n_func_basefuncs, n_qpoints)
    dNdX = fill(zero(Tensor{2,dim_s,T}) * T(NaN), n_func_basefuncs, n_qpoints)

    covar_base = fill(zero(Tensor{2,dim_s,T}) * T(NaN), n_qpoints)

    # Geometry interpolation
    n_geom_basefuncs = getnbasefunctions(geom_interpol)
    M    = fill(zero(T)          * T(NaN), n_geom_basefuncs, n_qpoints)
    dMdξ = fill(zero(Vec{dim_p,T}) * T(NaN), n_geom_basefuncs, n_qpoints)

    for (qp, ξ) in enumerate(quad_rule.points)
        basefunc_count = 1
        for basefunc in 1:getnbasefunctions(func_interpol)
            dNdξ_temp, N_temp = Tensors.gradient(ξ -> Ferrite.shape_value(func_interpol, ξ, basefunc), ξ, :all)
            for comp in 1:dim_s
                N_comp = zeros(T, dim_s)
                N_comp[comp] = N_temp

                N[basefunc_count, qp] = Vec{dim_s,T}((N_comp...,))

                dN_comp = zeros(T, dim_s, dim_s)

                dN_comp[comp, 1:dim_p] = dNdξ_temp
                dNdξ[basefunc_count, qp] = Tensor{2,dim_s,T}((dN_comp...,))
                basefunc_count += 1
            end
        end
        for i in 1:n_geom_basefuncs 
            dMdξ[i, qp], M[i, qp] = Tensors.gradient(ξ -> mid_surf_value(geom_interpol, ξ, i), ξ, :all)
        end
    end

    detJdA = fill(T(NaN), n_qpoints)
    R = fill(zero(Tensor{2,dim_s,T}) *T(NaN), n_qpoints)
    MM = Tensors.n_components(Tensors.get_base(eltype(R)))
    SurfaceVectorValues{dim_p,dim_s,T,MM,refshape,typeof(func_interpol)}(N, dNdξ, dNdX, R, detJdA, M, dMdξ, quad_rule, covar_base, func_interpol)
end

Ferrite.getngeobasefunctions(cv::SurfaceVectorValues{dim,dim_s}) where {dim,dim_s} = size(cv.N, 1) ÷ dim_s
Ferrite.getnbasefunctions(cv::SurfaceVectorValues) = size(cv.N, 1)

@inline getdetJdA(cv::SurfaceVectorValues, q_point::Int) = cv.detJdA[q_point]
@inline Ferrite.getdetJdV(cv::SurfaceVectorValues, q_point::Int) = getdetJdA(cv, q_point)

function Ferrite.reinit!(cv::SurfaceVectorValues{dim_p,dim_s}, x::AbstractVector{Vec{dim_s,T}}) where {dim_p,dim_s,T}
    n_geom_basefuncs = Ferrite.getngeobasefunctions(cv)
    @assert length(x) == n_geom_basefuncs


    @inbounds for i in 1:length(cv.qr.weights)
        w = cv.qr.weights[i]

        E = zeros(Vec{dim_s,T},dim_p)
        for j in 1:n_geom_basefuncs
            for d in 1:dim_p
                E[d] += cv.dMdξ[j,i][d] * x[j]
            end
        end
        D = Tensors.cross(E...) #in 2d cross-product is defined as cross(a::Vec{2}) = [-a[2], a[1]]
        
        #Rotation matrix and covariant vectors are similar becuase 
        _R = hcat((E./norm.(E))..., D/norm(D))
        _G = hcat(E...,D)

        cv.R[i] = Tensor{2,dim_s,T}(_R)
        cv.covar_base[i] = Tensor{2,dim_s,T}(_G)
        
        detJ = norm(D)#sqrt(det(cv.covar_base[i]))
        cv.detJdA[i] = detJ * w

        #Update dNdX
        #...not need at the moment

    end

end

function Ferrite.shape_value(fe_v::SurfaceVectorValues{dim,dim_s}, q_point::Int, i::Int) where {dim,T,dim_s}
    return fe_v.N[i,q_point]
end

function Ferrite.function_value(fe_v::SurfaceVectorValues{dim,dim_s}, q_point::Int, u::AbstractVector{T}, dof_range::AbstractVector{Int} = collect(1:length(u))) where {dim,T,dim_s}
    @assert length(u) == length(dof_range)
    val = zero(Vec{dim_s,T})
    for (i, j) in enumerate(dof_range)
        val += shape_value(fe_v, q_point, j) * u[i]
    end
    return val
end

