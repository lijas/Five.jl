Ferrite.cellid(f::Ferrite.BoundaryIndex) = f[1]

_ntupelcomponents(::Type{NTuple{N,Int}}) where N = N
Ferrite.nnodes(c::Type{<:Ferrite.AbstractCell}) = _ntupelcomponents(fieldtypes(c)[1])


_getquadraturerule(cv::CellValues) = cv.qr
_getinterpolation(cv::CellValues) = Ferrite.function_interpolation(cv)

function Ferrite.Dirichlet(;field::Symbol,set::Set{T},func::Function,dofs::Vector{Int}) where T
    return Ferrite.Dirichlet(field, set, func, dofs)
end

function getsubdofhandler(dh::DofHandler, cellid::Int)
    return dh.subdofhandlers[dh.cell_to_subdofhandler[cellid]]
end

Tensors.cross(v::Vec{2}) = Vec{2}((-v[2], v[1]))

##
#CellValues

function Ferrite.reference_coordinates(::Serendipity{dim,RefCube,0}) where dim
    return [zero(Vec{dim, Float64})]
end

function Ferrite.value(ip::Serendipity{dim,RefCube,0}, i::Int, ξ::Vec{dim}) where dim
    i == 1 && return 1
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

face_coordinate(dh::DofHandler{dim,T}, element::AbstractElement, faceindex::FaceIndex, location::Int) where {dim,T} = 
    _x_coordinate(dh, element, faceindex, location, Ferrite.faces)   
vertex_coordinate(dh::DofHandler{dim,T}, element::AbstractElement, vertexindex::FaceIndex, location::Int) where {dim,T} = 
    _x_coordinate(dh, element, vertexindex, location, Ferrite.vertices)

function _x_coordinate(dh::DofHandler{dim,T}, element, faceindex, location::Int, func::F) where {dim,T,F}
    xh = getcoordinates(dh.grid, faceindex[1])

    bcvalues = BCValues(get_field(element, :u).interpolation, default_interpolation(celltype(element)), func)
    #bcvalues = get_bcvalue(element, :u) 

    bcvalues.current_face[] = faceindex[2]

    x = spatial_coordinate(bcvalues,location, xh)

    return x
end

function faceset_to_nodes(grid, master_faceset)

    face_nodes = Int[]
    for faceindex in master_faceset
        cellid = faceindex[1]
        faceid = faceindex[2]

        facenodes = faces(getcells(grid, cellid))[faceid]
        push!(face_nodes, facenodes[1])
        push!(face_nodes, facenodes[2])
    end
    return unique(face_nodes)

end

function faceset_to_cellset(grid, master_faceset)

    face_nodes_set = Tuple{Int,Int}[]
    for faceindex in master_faceset
        cellid = faceindex[1]
        faceid = faceindex[2]

        facenodes = faces(getcells(grid, cellid))[faceid]
        push!(face_nodes_set, facenodes)
    end
    return face_nodes_set

end

function Base.getindex(v::AbstractArray, i::Vec{2,Int})
       return Vec{2,eltype(v)}((v[i[1]],v[i[2]]))
end

function Base.getindex(v::AbstractArray, i::Vec{3,Int})
       return Vec{3,eltype(v)}((v[i[1]],v[i[2]],v[i[3]]))
end

function second_derivative(ip::Interpolation{dim}, ξ::Vec{dim,T}) where {dim,T}
    [hessian(ξ -> value(ip, i, ξ), ξ) for i in 1:getnbasefunctions(ip)]
end

function generate_ooplane_quadraturerule(::Type{T}, zcoords::Vector{T}; nqp_per_layer::Int) where {T}
    
    nlayers = length(zcoords)-1

    points = Vector{Vec{1,T}}()
    weights = Vector{T}()
    oqr = QuadratureRule{1,RefCube}(nqp_per_layer)
    #copy the oo-plane integration to each layer
    addon = (last(zcoords) + first(zcoords))/2
    scale = (last(zcoords) - first(zcoords))/2
    zcoords = (zcoords.-addon)/scale
    for ilay in 1:nlayers
        addon = (zcoords[ilay+1] + zcoords[ilay])/2
        scale = (zcoords[ilay+1] - zcoords[ilay])/2
        for qp in 1:length(oqr.weights)
            new_z = oqr.points[qp]*scale .+ addon
            push!(points,  Vec(Tuple(new_z)))
            push!(weights, oqr.weights[qp]*scale)
        end
    end
    return QuadratureRule{1,RefCube,T}(weights,points)

end

function Ferrite.QuadratureRule{2,RefCube}(orders::NTuple{2,Int}) 
    T = Float64
    dim = 2
    
    qrs = [QuadratureRule{1,RefCube}(orders[d]) for d in 1:dim]

    newpoints = Vec{dim,T}[]
    newweights = T[]
        for i2 in 1:orders[2]
            for i1 in 1:orders[1]
                _v = (qrs[1].points[i1][1], qrs[2].points[i2][1])
                _w  = qrs[1].weights[i1] *  qrs[2].weights[i2] 
                push!(newpoints, Vec{dim,T}(_v))
                push!(newweights, _w)
            end
        end

    return QuadratureRule{dim,RefCube,T}(newweights, newpoints)
end

function Ferrite.QuadratureRule{3,RefCube}(orders::NTuple{3,Int}) 
    T = Float64
    dim = 3
    
    qrs = [QuadratureRule{1,RefCube}(orders[d]) for d in 1:dim]

    newpoints = Vec{dim,T}[]
    newweights = T[]
    for i3 in 1:orders[3]
        for i2 in 1:orders[2]
            for i1 in 1:orders[1]
                _v = (qrs[1].points[i1][1], qrs[2].points[i2][1], qrs[3].points[i3][1])
                _w  = qrs[1].weights[i1] *  qrs[2].weights[i2] * qrs[3].weights[i3]
                push!(newpoints, Vec{dim,T}(_v))
                push!(newweights, _w)
            end
        end
    end

    return QuadratureRule{dim,RefCube,T}(newweights, newpoints)
end

function merge_quadrules(qrs::Vararg{QuadratureRule,N}) where N
    @assert(N == 2 || N == 3)

    dim = N
    T = Float64
    #make this prettier:
    if N == 3
        newpoints = Vec{dim,T}[]
        newweights = T[]
        for i3 in 1:length(qrs[3].points)
            for i2 in 1:length(qrs[2].points)
                for i1 in 1:length(qrs[1].points)
                    _v = (qrs[1].points[i1][1], qrs[2].points[i2][1], qrs[3].points[i3][1])
                    _w  = qrs[1].weights[i1] *  qrs[2].weights[i2] * qrs[3].weights[i3]
                    push!(newpoints, Vec{dim,T}(_v))
                    push!(newweights, _w)
                end
            end
        end
    
        return QuadratureRule{dim,RefCube,T}(newweights, newpoints)
    elseif N == 2
        newpoints = Vec{dim,T}[]
        newweights = T[]
        for i2 in 1:length(qrs[2].points)
            for i1 in 1:length(qrs[1].points)
                _v = (qrs[1].points[i1][1], qrs[2].points[i2][1])
                _w  = qrs[1].weights[i1] *  qrs[2].weights[i2] 
                push!(newpoints, Vec{dim,T}(_v))
                push!(newweights, _w)
            end
        end
    
        return QuadratureRule{dim,RefCube,T}(newweights, newpoints)
    end

end

#TODO: Clean this
##Adding edge set to dbc
##
##
#
function add_edge!(ch::ConstraintHandler, dbc::Dirichlet)
    Ferrite.dbc_check(ch, dbc)
    field_idx = Ferrite.find_field(ch.dh, dbc.field_name)
    # Extract stuff for the field
    interpolation =Ferrite.getfieldinterpolation(ch.dh, field_idx)#ch.dh.field_interpolations[field_idx]
    field_dim =Ferrite. getfielddim(ch.dh, field_idx)#ch.dh.field_dims[field_idx] # TODO: I think we don't need to extract these here ...
    bcvalue = Ferrite.getbcvalue(ch.dh, field_idx)
    _add_edge!(ch, dbc, dbc.faces, interpolation, field_dim, Ferrite.field_offset(ch.dh, dbc.field_name), bcvalue)
    return ch
end

function _add_edge!(ch::ConstraintHandler, dbc::Dirichlet, bcfaces::Set{Tuple{Int,Int}}, interpolation::Interpolation, field_dim::Int, offset::Int, bcvalue::Ferrite.BCValues)
    # calculate which local dof index live on each face
    # face `i` have dofs `local_face_dofs[local_face_dofs_offset[i]:local_face_dofs_offset[i+1]-1]
    local_face_dofs = Int[]
    local_face_dofs_offset = Int[1]
    for (i, face) in enumerate(Ferrite.edges(interpolation))
        for fdof in face, d in 1:field_dim
            if d ∈ dbc.components # skip unless this component should be constrained
                push!(local_face_dofs, (fdof-1)*field_dim + d + offset)
            end
        end
        push!(local_face_dofs_offset, length(local_face_dofs) + 1)
    end
    Ferrite.copy!!(dbc.local_face_dofs, local_face_dofs)
    Ferrite.copy!!(dbc.local_face_dofs_offset, local_face_dofs_offset)

    # loop over all the faces in the set and add the global dofs to `constrained_dofs`
    constrained_dofs = Int[]
    #_celldofs = fill(0, ndofs_per_cell(ch.dh))
    for (cellidx, faceidx) in bcfaces
        _celldofs = fill(0, ndofs_per_cell(ch.dh, cellidx))
        celldofs!(_celldofs, ch.dh, cellidx) # extract the dofs for this cell
        r = local_face_dofs_offset[faceidx]:(local_face_dofs_offset[faceidx+1]-1)
        append!(constrained_dofs, _celldofs[local_face_dofs[r]]) # TODO: for-loop over r and simply push! to ch.prescribed_dofs
        @debug println("adding dofs $(_celldofs[local_face_dofs[r]]) to dbc")
    end

    # save it to the ConstraintHandler
    push!(ch.dbcs, dbc)
    push!(ch.bcvalues, bcvalue)
    append!(ch.prescribed_dofs, constrained_dofs)
end

#=function Ferrite.spatial_coordinate(bcv::Ferrite.BCValues, q_point::Int, xh::AbstractVector{Vec{dim,T}}) where {dim,T}
    return zero(Vec{dim,T})
end=#


#Overwrite Ferrite.function_value for scalar cellvalues to interpolate any type
#=
function Ferrite.function_value(fe_v::Ferrite.ScalarValues{dim}, q_point::Int, u::AbstractVector{T}, dof_range::AbstractVector{Int} = collect(1:length(u))) where {dim,T}
    n_base_funcs = Ferrite.getn_scalarbasefunctions(fe_v)
    #isa(fe_v, VectorValues) && (n_base_funcs *= dim)
    @assert length(dof_range) == n_base_funcs
    @boundscheck checkbounds(u, dof_range)
    val = zero(T)#val = zero(_valuetype(fe_v, u))
    @inbounds for (i, j) in enumerate(dof_range)
        val += shape_value(fe_v, q_point, i) * u[j]
    end
    return val
end

function function_geometric_derivative(fe_v::Ferrite.ScalarValues{dim}, q_point::Int, u::AbstractVector{T}, d::Int, dof_range::AbstractVector{Int} = collect(1:length(u))) where {dim,T}
    n_base_funcs = Ferrite.getn_scalarbasefunctions(fe_v)
    @assert length(dof_range) == n_base_funcs
    @boundscheck checkbounds(u, dof_range)
    grad = zero(T)
    @inbounds for (i, j) in enumerate(dof_range)
        grad += fe_v.dMdξ[i,q_point][d] * u[j]
    end
    return grad
end

function function_shape_derivative(fe_v::Ferrite.ScalarValues{dim_p}, q_point::Int, x::Vector{Vec{dim_s,T}}, u::AbstractVector{T2}, a::Vec{dim_s,T}, b::Vec{dim_s,T}, dof_order::AbstractVector{Int} = collect(1:length(u))) where {dim_p,dim_s,T,T2}
    n_base_funcs = Ferrite.getn_scalarbasefunctions(fe_v)
    @assert length(dof_order) == n_base_funcs
    @boundscheck checkbounds(u, dof_order)

    dudξ = zeros(T2,dim_p)
    dxdξ = zeros(Vec{dim_s,T},dim_p)

    @inbounds for (i, j) in enumerate(dof_order)
        for d in 1:dim_p
            dxdξ[d] += fe_v.dMdξ[i,q_point][d] * x[i]
            dudξ[d] += fe_v.dMdξ[i,q_point][d] * u[j]
        end
    end

    #ab = [a,b]
    #A = [ab[i]⋅ab[j] for i in 1:dim_p, j in 1:dim_p]
    #B = [ab[i]⋅dxdξ[j] for i in 1:dim_p, j in 1:dim_p]

    #dudξ̂ ::Vector{T2} = A'*inv(B') * dudξ

    return dudξ#dudξ̂ 

end


function Ferrite.spatial_coordinate(fe_v::Ferrite.ScalarValues{dim_p}, q_point::Int, x::AbstractVector{Vec{dim_s,T}}) where {dim_p,dim_s,T}
    n_base_funcs = Ferrite.getngeobasefunctions(fe_v)
    @assert length(x) == n_base_funcs
    vec = zero(Vec{dim_s,T})
    @inbounds for i in 1:n_base_funcs
        vec += Ferrite.geometric_value(fe_v, q_point, i) * x[i]
    end
    return vec
end
##
=#

function gridmerge(grids::Vararg{Grid{dim,C,T}}) where dim where T where C
    
    nodes_new = Node{dim,T}[]
    cells_new = Ferrite.AbstractCell[]

    faceset_new = Dict{String, Set{FaceIndex}}()
    vertexset_new = Dict{String, Set{VertexIndex}}()
    edgeset_new = Dict{String, Set{EdgeIndex}}()
    nodeset_new = Dict{String, Set{Int}}()
    cellset_new = Dict{String, Set{Int}}()

    for (igrid, grid) in enumerate(grids)
        nodeoffset = length(nodes_new)
        celloffset = length(cells_new)

        #For all cells in grid2, increment the nodeids
        for cell in grid.cells
            N = length(cell.nodes)
            M = nfaces(cell)
            offset_nodes = [n + nodeoffset for n in cell.nodes]
            offset_cell  = Cell{dim,N,M}(NTuple{N,Int}(tuple(offset_nodes...)));
            push!(cells_new, offset_cell)
        end
        append!(nodes_new, grid.nodes)

        #Update sets
        function updata_set(indexset, new_indexset, INDEX)
            grid2_new_facesets = Dict{String, Set{INDEX}}()
            for (facesetname, faceset) in getproperty(grid, indexset)
                tmp_faceset = Set{INDEX}()
                for (cellid, faceidx) in faceset
                    push!(tmp_faceset, INDEX(cellid + celloffset, faceidx))
                end
                new_indexset[facesetname * string(igrid)] =  tmp_faceset
            end
        end

        updata_set(:facesets, faceset_new, FaceIndex)
        updata_set(:vertexsets, vertexset_new, VertexIndex)
        updata_set(:edgesets, edgeset_new, EdgeIndex)

        for (nodesetname, nodeset) in grid.nodesets
            tmp_nodeset = Set{Int}()
            for nodeid in nodeset
                push!(tmp_nodeset, nodeid + nodeoffset)
            end
            nodeset_new[nodesetname * String(igrid)] =  tmp_nodeset
        end

        for (setname, cellset) in grid.cellsets
            tmp_nodeset = Set{Int}()
            for cellid in cellset
                push!(tmp_nodeset, cellid + celloffset)
            end
            cellset_new[setname * String(igrid)] =  tmp_nodeset
        end

    end
    bm = spzeros(Bool, 0, 0)
    
    return Grid(cells_new, nodes_new, cellset_new, nodeset_new, faceset_new, edgeset_new, vertexset_new, bm)

end

function WriteVTK.vtk_cell_data(
    vtk::WriteVTK.DatasetFile,
    data::Vector{S},
    name::AbstractString
    ) where {O, D, T, M, S <: Union{Tensor{O, D, T, M}, SymmetricTensor{O, D, T, M}}}
    ncells = length(data)
    out = zeros(T, M, ncells)
    for i in 1:ncells
        Ferrite.toparaview!(@view(out[:, i]), data[i])
    end
    return vtk_cell_data(vtk, out, name; component_names=component_names(S))
end

function component_names(::Type{S}) where S
    names =
        S <:             Vec{1}   ? ["x"] :
        S <:             Vec      ? ["x", "y", "z"] : # Pad 2D Vec to 3D
        S <:          Tensor{2,1} ? ["xx"] :
        S <: SymmetricTensor{2,1} ? ["xx"] :
        S <:          Tensor{2,2} ? ["xx", "yy", "xy", "yx"] :
        S <: SymmetricTensor{2,2} ? ["xx", "yy", "xy"] :
        S <:          Tensor{2,3} ? ["xx", "yy", "zz", "yz", "xz", "xy", "zy", "zx", "yx"] :
        S <: SymmetricTensor{2,3} ? ["xx", "yy", "zz", "yz", "xz", "xy"] :
                                    nothing
    return names
end


export disassemble!
Base.@propagate_inbounds function disassemble!(ue::AbstractVector, u::AbstractVector, dofs::AbstractVector{Int})
    Base.@boundscheck checkbounds(u, dofs)
    # @inbounds for i in eachindex(ue, dofs) # Slow on Julia 1.6 (JuliaLang/julia#40267)
        Base.@inbounds for i in eachindex(ue)
        ue[i] = u[dofs[i]]
    end
    return ue
end

function Ferrite._check_same_celltype(grid::Ferrite.AbstractGrid, boundaryset::AbstractVector{<:Ferrite.BoundaryIndex})
    cellid, faceid = first(boundaryset)
    celltype = typeof(grid.cells[cellid])
    for (cellid,bid) in boundaryset
        if celltype != typeof(grid.cells[cellid])
            error("The cells in your boundary-set are not all of the same celltype.")
        end
    end
end
