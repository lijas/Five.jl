export CohesiveCell

"""

"""

struct CohesiveZoneInterpolation{dim,order,I<:Interpolation} <: Ferrite.Interpolation{dim,RefCube,order} 
    interpolation::I
    function CohesiveZoneInterpolation(interpolation)
        dim = Ferrite.getdim(interpolation) + 1
        order = Ferrite.getorder(interpolation)
        new{dim,order,typeof(interpolation)}(interpolation)
    end
end

Ferrite.getnbasefunctions(ip::CohesiveZoneInterpolation) = getnbasefunctions(ip.interpolation)*2
Ferrite.nvertexdofs(::CohesiveZoneInterpolation) = 1
Ferrite.nfacedofs(ip::CohesiveZoneInterpolation) = Ferrite.ncelldofs(ip.interpolation)

_mapper(::Lagrange{1,RefCube,1}, _i::Int) = [1,2,1,2][_i]
_mapper(::Lagrange{1,RefCube,2}, _i::Int) = [1,2,1,2,3,3][_i]
_mapper(::Lagrange{2,RefCube,1}, _i::Int) = [1,2,3,4,1,2,3,4][_i]
function mid_surf_value(ip::CohesiveZoneInterpolation, _i::Int, ξ::Vec{dim_p}) where dim_p
    _i<=getnbasefunctions(ip) || throw(ArgumentError("no shape function $_i for interpolation $ip"))
    i = _mapper(ip.interpolation, _i)

    return Ferrite.value(ip.interpolation, i, ξ)*0.5
end

_sign_mapper(::Lagrange{1,RefCube,1}, _i::Int) = [-1,-1,+1,+1][_i]
_sign_mapper(::Lagrange{1,RefCube,2}, _i::Int) = [-1,-1,+1,+1,-1,+1][_i]
_sign_mapper(::Lagrange{2,RefCube,1}, _i::Int) = [-1,-1,-1,-1,+1,+1,+1,+1][_i]
function Ferrite.value(ip::CohesiveZoneInterpolation, _i::Int, ξ::Vec{dim_p}) where dim_p
    _i<=getnbasefunctions(ip) || throw(ArgumentError("no shape function $i for interpolation $ip"))
    sign = _sign_mapper(ip.interpolation, _i)
    i = _mapper(ip.interpolation, _i)
    return sign*Ferrite.value(ip.interpolation, i, ξ)
end

Ferrite._mass_qr(::CohesiveZoneInterpolation{2,1,Lagrange{1,RefCube,1}}) = QuadratureRule{1,RefCube}(1)

"""

"""

struct CohesiveCell{dim,N,M} <: Ferrite.AbstractCell{dim,N,M}
    nodes::NTuple{N,Int}
end

Ferrite.vertices(c::CohesiveCell{2,N,2}) where N = (c.nodes[1], c.nodes[2], c.nodes[3], c.nodes[4])
Ferrite.faces(c::CohesiveCell{2,N,2}) where N    = ((c.nodes[1], c.nodes[2]), (c.nodes[3], c.nodes[4])) 

Ferrite.vertices(c::CohesiveCell{3,8,2}) = (c.nodes[1], c.nodes[2], c.nodes[3], c.nodes[4], c.nodes[5], c.nodes[6], c.nodes[7], c.nodes[8])
Ferrite.faces(c::CohesiveCell{3,8,2}) = ((c.nodes[1], c.nodes[2], c.nodes[3], c.nodes[4]), (c.nodes[5], c.nodes[6], c.nodes[7], c.nodes[8])) 

Ferrite.default_interpolation(::Type{CohesiveCell{2,4,2}}) = CohesiveZoneInterpolation(Lagrange{1,RefCube,1}())
Ferrite.default_interpolation(::Type{CohesiveCell{2,6,2}}) = CohesiveZoneInterpolation(Lagrange{1,RefCube,2}())
Ferrite.default_interpolation(::Type{CohesiveCell{3,8,2}}) = CohesiveZoneInterpolation(Lagrange{2,RefCube,1}())

Ferrite.cell_to_vtkcell(::Type{CohesiveCell{2,4,2}}) = Ferrite.VTKCellTypes.VTK_QUAD
Ferrite.cell_to_vtkcell(::Type{CohesiveCell{3,8,2}}) = Ferrite.VTKCellTypes.VTK_HEXAHEDRON

#Ferrite.nodes_to_vtkorder(cell::CohesiveCell{2,4,2}) = cell.nodes[[1,2,4,3]]

#Ferrite.cell_to_vtkcell(::Type{CohesiveCell{2,4,2}}) = Ferrite.VTKCellTypes.VTK_LINE
#Ferrite.cell_to_vtkcell(::Type{CohesiveCell{3,8,2}}) = Ferrite.VTKCellTypes.VTK_QUAD

"""
    SurfaceVectorValues

Contains shape values for cohesive zone elements
"""
struct SurfaceVectorValues{dim_p,dim_s,T<:Real,M2,refshape<:Ferrite.AbstractRefShape} <: Ferrite.CellValues{dim_s,T,refshape}
    N::Matrix{Vec{dim_s,T}}
    dNdξ::Matrix{Tensor{2,dim_s,T,M2}}
    dNdx::Matrix{Tensor{2,dim_s,T,M2}}
    R::Vector{Tensor{2,dim_s,T,M2}}
    detJdA::Vector{T}
    M::Matrix{T}  # Shape values for geometric interp
    dMdξ::Matrix{Vec{dim_p,T}}
    qr::QuadratureRule{dim_p,refshape,T}
    covar_base::Vector{Tensor{2,dim_s,T}}
end

Ferrite.getnquadpoints(cv::SurfaceVectorValues) = length(cv.qr.weights)

function SurfaceVectorValues(quad_rule::QuadratureRule, func_interpol::CohesiveZoneInterpolation, geom_interpol::CohesiveZoneInterpolation=func_interpol)
    SurfaceVectorValues(Float64, quad_rule, func_interpol, geom_interpol)
end

@inline getR(cv::SurfaceVectorValues, qp::Int) = cv.R[qp]

function SurfaceVectorValues(::Type{T},
    quad_rule::QuadratureRule{dim_p,RefCube},
    func_interpol::CohesiveZoneInterpolation{dim_s},
    geom_interpol::CohesiveZoneInterpolation{dim_s}=func_interpol) where {dim_p,dim_s,T}

    @assert Ferrite.getdim(func_interpol) == Ferrite.getdim(geom_interpol)
    @assert Ferrite.getrefshape(func_interpol) == Ferrite.getrefshape(geom_interpol) == RefCube
    n_qpoints = length(getweights(quad_rule))

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
            dNdξ_temp, N_temp = Tensors.gradient(ξ -> Ferrite.value(func_interpol, basefunc, ξ), ξ, :all)
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
            dMdξ[i, qp], M[i, qp] = Tensors.gradient(ξ -> mid_surf_value(geom_interpol, i, ξ), ξ, :all)
        end
    end

    detJdA = fill(T(NaN), n_qpoints)
    R = fill(zero(Tensor{2,dim_s,T}) *T(NaN), n_qpoints)
    MM = Tensors.n_components(Tensors.get_base(eltype(R)))
    SurfaceVectorValues{dim_p,dim_s,T,MM,RefCube}(N, dNdξ, dNdX, R, detJdA, M, dMdξ, quad_rule, covar_base)
end

Ferrite.getn_scalarbasefunctions(cv::SurfaceVectorValues{dim,dim_s}) where {dim,dim_s} = size(cv.N, 1) ÷ dim_s

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
        D = cross(E...) #in 2d cross-product is defined as cross(a::Vec{2}) = [-a[2], a[1]]
        
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


function Ferrite.function_value(fe_v::SurfaceVectorValues{dim,dim_s}, q_point::Int, u::AbstractVector{T}, dof_range::AbstractVector{Int} = collect(1:length(u))) where {dim,T,dim_s}
    n_base_funcs = Ferrite.getn_scalarbasefunctions(fe_v)
    n_base_funcs *= dim_s
    
    @assert length(u) == length(dof_range)
    val = zero(Vec{dim_s,T})
    for (i, j) in enumerate(dof_range)
        val += shape_value(fe_v, q_point, j) * u[i]
    end
    return val
end

