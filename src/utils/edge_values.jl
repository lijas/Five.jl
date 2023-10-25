export EdgeVectorValues, EdgeValues

struct EdgeVectorValues{dim,T<:Real,refshape<:Ferrite.AbstractRefShape,M} 
    N::Array{Vec{dim,T},3}
    dNdx::Array{Tensor{2,dim,T,M},3}
    dNdξ::Array{Tensor{2,dim,T,M},3}
    detJdV::Matrix{T}
    normals::Vector{Vec{dim,T}}
    M::Array{T,3}
    dMdξ::Array{Vec{dim,T},3}
    qr_weights::Vector{T}
    current_face::Ferrite.ScalarWrapper{Int}
end

function EdgeVectorValues(quad_rule::QuadratureRule, func_interpol::Ferrite.Interpolation, geom_interpol::Ferrite.Interpolation=func_interpol)
    EdgeVectorValues(Float64, quad_rule, func_interpol, geom_interpol)
end

function EdgeVectorValues(::Type{T}, quad_rule::QuadratureRule{dim_qr,shape}, func_interpol::Ferrite.Interpolation,
        geom_interpol::Ferrite.Interpolation=func_interpol) where {dim_qr,T,shape<:Ferrite.Ferrite.AbstractRefShape}

    @assert Ferrite.getdim(func_interpol) == Ferrite.getdim(geom_interpol)
    @assert Ferrite.getrefshape(func_interpol) == Ferrite.getrefshape(geom_interpol) == shape
    n_qpoints = length(getweights(quad_rule))
    dim = dim_qr + 2

    edge_quad_rule = create_edge_quad_rule(quad_rule, func_interpol)
    n_faces = length(edge_quad_rule)

    # Normals
    normals = zeros(Vec{dim,T}, n_qpoints)

    # Function interpolation
    n_func_basefuncs = getnbasefunctions(func_interpol) * dim
    N    = fill(zero(Vec{dim,T})      * T(NaN), n_func_basefuncs, n_qpoints, n_faces)
    dNdx = fill(zero(Tensor{2,dim,T}) * T(NaN), n_func_basefuncs, n_qpoints, n_faces)
    dNdξ = fill(zero(Tensor{2,dim,T}) * T(NaN), n_func_basefuncs, n_qpoints, n_faces)

    # Geometry interpolation
    n_geom_basefuncs = getnbasefunctions(geom_interpol)
    M    = fill(zero(T)          * T(NaN), n_geom_basefuncs, n_qpoints, n_faces)
    dMdξ = fill(zero(Vec{dim,T}) * T(NaN), n_geom_basefuncs, n_qpoints, n_faces)

    for face in 1:n_faces, (qp, ξ) in enumerate(edge_quad_rule[face].points)
        basefunc_count = 1
        for basefunc in 1:getnbasefunctions(func_interpol)
            dNdξ_temp, N_temp = gradient(ξ -> Ferrite.value(func_interpol, basefunc, ξ), ξ, :all)
            for comp in 1:dim
                N_comp = zeros(T, dim)
                N_comp[comp] = N_temp
                N[basefunc_count, qp, face] = Vec{dim,T}((N_comp...,))

                dN_comp = zeros(T, dim, dim)
                dN_comp[comp, :] = dNdξ_temp
                dNdξ[basefunc_count, qp, face] = Tensor{2,dim,T}((dN_comp...,))
                basefunc_count += 1
            end
        end
        for basefunc in 1:n_geom_basefuncs
            dMdξ[basefunc, qp, face], M[basefunc, qp, face] = gradient(ξ -> Ferrite.value(geom_interpol, basefunc, ξ), ξ, :all)
        end
    end

    detJdV = fill(T(NaN), n_qpoints, n_faces)
    MM = Tensors.n_components(Tensors.get_base(eltype(dNdx)))

    EdgeVectorValues{dim,T,shape,MM}(N, dNdx, dNdξ, detJdV, normals, M, dMdξ, quad_rule.weights, Ferrite.ScalarWrapper(0))
end

function Ferrite.reinit!(fv::EdgeValues{dim}, x::AbstractVector{Vec{dim,T}}, face::Int) where {dim,T}
    n_geom_basefuncs = Ferrite.getngeobasefunctions(fv)
    n_func_basefuncs = Ferrite.getn_scalarbasefunctions(fv)
    @assert length(x) == n_geom_basefuncs
    isa(fv, EdgeVectorValues) && (n_func_basefuncs *= dim)

    fv.current_face[] = face
    cb = getcurrentedge(fv)

    @inbounds for i in 1:length(fv.qr_weights)
        w = fv.qr_weights[i]
        fefv_J = zero(Tensor{2,dim})
        for j in 1:n_geom_basefuncs
            fefv_J += x[j] ⊗ fv.dMdξ[j, i, cb]
        end
        weight_norm = weighted_normal(fefv_J, fv, cb)
        fv.normals[i] = weight_norm / norm(weight_norm)
        detJ = norm(weight_norm)
        
        detJ > 0.0 || throw(ArgumentError("det(J) is not positive: det(J) = $(detJ)"))
        fv.detJdV[i, cb] = detJ * w
        Jinv = inv(fefv_J)
        for j in 1:n_func_basefuncs
            fv.dNdx[j, i, cb] = fv.dNdξ[j, i, cb] ⋅ Jinv
        end
    end
end

getcurrentedge(fv::EdgeValues) = fv.current_face[]
Ferrite.getnormal(fv::EdgeValues, qp::Int) = fv.normals[qp]
Ferrite.getdetJdV(bv::EdgeValues, q_point::Int) = bv.detJdV[q_point, bv.current_face[]]

Ferrite.getn_scalarbasefunctions(cv::EdgeVectorValues{dim}) where {dim} = size(cv.N, 1) ÷ dim

##################
# All RefCube 3D #
##################
function create_edge_quad_rule(quad_rule::QuadratureRule{1,shape,T}, ::Ferrite.Interpolation{3,shape}) where {T,shape<:RefCube}
    w = getweights(quad_rule)
    p = getpoints(quad_rule)
    n_points = length(w)
    edge_quad_rule = QuadratureRule{3,shape,T}[]

    # Bottom
    # Edge 1
    new_points = [Vec{3,T}((p[i][1], -one(T), -one(T))) for i in 1:n_points] # ξ = t, η = -1, ζ = -1
    push!(edge_quad_rule, QuadratureRule{3,shape,T}(w, new_points))
    # Edge 2
    new_points = [Vec{3,T}((one(T), p[i][1], -one(T))) for i in 1:n_points] # ξ = 1, η = s, ζ = -1
    push!(edge_quad_rule, QuadratureRule{3,shape,T}(w, new_points))
    # Edge 3
    new_points = [Vec{3,T}((p[i][1], one(T), -one(T))) for i in 1:n_points] # ξ = t, η = s, ζ = -1
    push!(edge_quad_rule, QuadratureRule{3,shape,T}(w, new_points))
    # Edge 4
    new_points = [Vec{3,T}((-one(T), p[i][1], -one(T))) for i in 1:n_points] # ξ = t, η = s, ζ = -1
    push!(edge_quad_rule, QuadratureRule{3,shape,T}(w, new_points))            
    
    # Top
    # Edge 1
    new_points = [Vec{3,T}((p[i][1], -one(T), one(T))) for i in 1:n_points] # ξ = t, η = s, ζ = -1
    push!(edge_quad_rule, QuadratureRule{3,shape,T}(w, new_points))
    # Edge 2
    new_points = [Vec{3,T}((one(T), p[i][1], one(T))) for i in 1:n_points] # ξ = t, η = s, ζ = -1
    push!(edge_quad_rule, QuadratureRule{3,shape,T}(w, new_points))
    # Edge 3
    new_points = [Vec{3,T}((p[i][1], one(T), one(T))) for i in 1:n_points] # ξ = t, η = s, ζ = -1
    push!(edge_quad_rule, QuadratureRule{3,shape,T}(w, new_points))
    # Edge 4
    new_points = [Vec{3,T}((-one(T), p[i][1], one(T))) for i in 1:n_points] # ξ = t, η = s, ζ = -1
    push!(edge_quad_rule, QuadratureRule{3,shape,T}(w, new_points))     

    # Vertical edges
    # Edge 1
    new_points = [Vec{3,T}((-one(T), -one(T), p[i][1])) for i in 1:n_points] # ξ = t, η = s, ζ = -1
    push!(edge_quad_rule, QuadratureRule{3,shape,T}(w, new_points))
    # Edge 2
    new_points = [Vec{3,T}((one(T), -one(T), p[i][1])) for i in 1:n_points] # ξ = t, η = s, ζ = -1
    push!(edge_quad_rule, QuadratureRule{3,shape,T}(w, new_points))
    # Edge 3
    new_points = [Vec{3,T}((one(T), one(T), p[i][1])) for i in 1:n_points] # ξ = t, η = s, ζ = -1
    push!(edge_quad_rule, QuadratureRule{3,shape,T}(w, new_points))
    # Edge 4
    new_points = [Vec{3,T}((-one(T), one(T), p[i][1])) for i in 1:n_points] # ξ = t, η = s, ζ = -1
    push!(edge_quad_rule, QuadratureRule{3,shape,T}(w, new_points)) 

    return edge_quad_rule
end

function weighted_normal(J::Tensor{2,3}, ::EdgeValues{3,T,RefCube}, edge::Int) where {T}
    @inbounds begin
        edge == 1 && return J[:,1]
        edge == 2 && return J[:,2]
        edge == 3 && return J[:,1]
        edge == 4 && return J[:,2]
        edge == 5 && return J[:,1]
        edge == 6 && return J[:,2]
        edge == 7 && return J[:,1]
        edge == 8 && return J[:,2]
        edge == 9 && return J[:,3]
        edge == 10 && return J[:,3]
        edge == 11 && return J[:,3]
        edge == 12 && return J[:,3]
    end
    throw(ArgumentError("unknown edge number: $edge"))
end

@inline Ferrite.shape_value(bv::EdgeValues, q_point::Int, base_func::Int) = bv.N[base_func, q_point, bv.current_face[]]
Base.@pure Ferrite._valuetype(::EdgeValues{dim}, ::AbstractVector{T}) where {dim,T} = Vec{dim,T}

function Ferrite.function_value(fe_v::EdgeValues{dim}, q_point::Int, u::AbstractVector{T}, dof_range::UnitRange = 1:length(u)) where {dim,T}
    n_base_funcs = Ferrite.getn_scalarbasefunctions(fe_v)
    isa(fe_v, EdgeVectorValues) && (n_base_funcs *= dim)
    @assert length(dof_range) == n_base_funcs
    @boundscheck checkbounds(u, dof_range)
    val = zero(Ferrite._valuetype(fe_v, u))
    @inbounds for (i, j) in enumerate(dof_range)
        val += shape_value(fe_v, q_point, i) * u[j]
    end
    return val
end

#=
function testedge()
    N = 10
    L = 2.0
    left = zero(Vec{3})
    right = Vec{3}((10,2.0,3.0))
    grid = generate_grid(Hexahedron, (N, N, N), left, right)

    addedgeset!(grid, "edge1", (x)->x[1]≈0.0 && x[3]≈0.0)
    addedgeset!(grid, "edge2", (x)->x[2]≈0.0 && x[3]≈0.0)
    addedgeset!(grid, "edge3", (x)->x[1]≈0.0 && x[2]≈0.0)

    ip = Lagrange{3,RefCube,1}()

    qr = QuadratureRule{1,RefCube}(2)
    ev = EdgeVectorValues(qr, ip)

    V=0
    for (cellid, edgeidx) in getedgeset(grid, "edge1")
        x = getcoordinates(grid, cellid)
        reinit!(ev, x, edgeidx)

        for qp in 1:getnquadpoints(ev)
            V += getdetJdV(ev, qp)
        end
    end
    @show V

    V=0
    for (cellid, edgeidx) in getedgeset(grid, "edge2")
        x = getcoordinates(grid, cellid)
        reinit!(ev, x, edgeidx)

        for qp in 1:getnquadpoints(ev)
            V += getdetJdV(ev, qp)
        end
    end
    @show V

    V=0
    for (cellid, edgeidx) in getedgeset(grid, "edge3")
        x = getcoordinates(grid, cellid)
        reinit!(ev, x, edgeidx)

        for qp in 1:getnquadpoints(ev)
            V += getdetJdV(ev, qp)
        end
    end
    @show V
end

=#