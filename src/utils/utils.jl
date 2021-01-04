export nodeset_to_vertexset


"""
nodeset_to_vertexset(grid::AbstractGrid, nodeset::Set{Int}) 

Given a nodeset, convert it to a vertexset
"""
function nodeset_to_vertexset(grid::JuAFEM.AbstractGrid, nodeset::Set{Int}) 
    vertexset = Set{VertexIndex}()
    for nodeid in nodeset
        for (cellid, cell) in enumerate(grid.cells)
            for (i, nodeid2) in enumerate(cell.nodes)
                if nodeid == nodeid2
                    push!(vertexset, VertexIndex(cellid,i))
                    break
                end    
            end
        end
    end
    return vertexset
end


function indexdofs(dh, element::FieldHandler, index::Index, field_name::Symbol, components::Vector{Int})
    _dofs_on_something(dh, element, index, field_name, components, JuAFEM.getgeometryfunction(typeof(index)))
end
dofs_on_edge(dh, element::FieldHandler, faceindex::FaceIndex, field_name::Symbol, components::Vector{Int}) = 
    _dofs_on_something(dh, element, faceindex, field_name, components, JuAFEM.edges)
dofs_on_face(dh, element::FieldHandler, faceindex::FaceIndex, field_name::Symbol, components::Vector{Int}) = 
    _dofs_on_something(dh, element, faceindex, field_name, components, JuAFEM.faces)
dofs_on_vertex(dh, element::FieldHandler, vertexindex::VertexIndex, field_name::Symbol, components::Vector{Int}) = 
    _dofs_on_something(dh, element, vertexindex, field_name, components, JuAFEM.vertices)
function _dofs_on_something(dh::JuAFEM.AbstractDofHandler, fh::FieldHandler, faceindex, field_name::Symbol, components::Vector{Int}, faces_or_vertices::F) where {dim,T,F}
    
    #components = 1:dim
    #field =  element.fields[1]
    fieldidx = JuAFEM.find_field(fh, field_name)
    field = fh.fields[fieldidx]
    offset = JuAFEM.field_offset(fh, field_name)
    local_face_dofs = Int[]
    face = faces_or_vertices(field.interpolation)[faceindex[2]]
    
    for fdof in face, d in 1:field.dim
        if d ∈ components # skip unless this component should be constrained
            push!(local_face_dofs, (fdof-1)*field.dim + d + offset)
        end
    end
    
    # loop over all the faces in the set and add the global dofs to `constrained_dofs`
    _celldofs = celldofs(dh, faceindex[1])
    
    return _celldofs[local_face_dofs] # TODO: for-loop over r and simply push! to ch.prescribed_dofs

end

function get_node_dofs(dh::DofHandler{dim}, field_name = :u, ncomps = 2) where dim

    nnodes = getnnodes(dh.grid)
    node_dofs = zeros(Int, ncomps, nnodes)
    visited = falses(nnodes)
    for (ie, element) in enumerate(dh.elements)

        field = get_field(element, field_name)
        if field == nothing 
            continue
        end

        for cellidx in dh.elementcells[ie]
            cell = getcells(dh.grid, cellidx)

            field_dim = field.dim
            interpol_points = getnbasefunctions(field.interpolation)
            _celldofs = fill(0, ndofs(element))
            offset = field_offset(element, field_name)

            celldofs!(_celldofs, dh, cellidx) # update the dofs for this cell
            for idx in 1:length(cell.nodes)#min(interpol_points, length(cell.nodes))
                node = cell.nodes[idx]
                if !visited[node]
                    noderange = (offset + (idx-1)*field_dim + 1):(offset + idx*field_dim) # the dofs in this node
                    for (i,c) in enumerate(1:ncomps)
                        node_dofs[i,node] = _celldofs[noderange[c]]
                    end
                    visited[node] = true
                end
            end
        end
    end
    return node_dofs
end

#Gets all initial coordinates (:u and :xyθ)
function get_x0(dh::JuAFEM.AbstractDofHandler)

    T = getT(dh)
    dim = JuAFEM.getdim(dh)

    nnodes = getnnodes(dh.grid)
    x0 = zeros(T, ndofs(dh))
    visited = falses(nnodes)
    for (ie, element) in enumerate(dh.fieldhandlers)

        idx = findfirst(i->i == :u, JuAFEM.getfieldnames(element))
        field = element.fields[idx]
        if field === nothing 
            continue
            #field = get_field(element, :xyθ)
            #if field == nothing 
            #    continue
            #end
        end
        field_dim = field.dim
        offset = JuAFEM.field_offset(element, field.name)
        _celldofs = Int[]

        for cellidx in element.cellset
            resize!(_celldofs, ndofs_per_cell(dh, cellidx))
            cell = getcells(dh.grid, cellidx)
            celldofs!(_celldofs, dh, cellidx) # update the dofs for this cell
            for idx in 1:length(cell.nodes)#min(interpol_points, length(cell.nodes))
                node = cell.nodes[idx]
                if !visited[node]

                    noderange = (offset + (idx-1)*field_dim + 1):(offset + idx*field_dim) # the dofs in this node

                    for c = 1:dim
                        d = _celldofs[noderange[c]]
                        x0[d] = dh.grid.nodes[node].x[c]
                    end
                    visited[node] = true
                end
            end
        end
    end
    return x0

end


function numdiff(f!::Function, ue_interface::AbstractVector{T}, nms::AbstractVector{<:AbstractMaterialState}, ms::AbstractVector{<:AbstractMaterialState}) where T

    h = 1e-7
    ndofs = length(ue_interface)
    ke_temp = zeros(T, ndofs, ndofs)
    fdu = zeros(T, ndofs)
    fu = zeros(T, ndofs)
    f₊ = zeros(T, ndofs)
    f₋ = zeros(T, ndofs)
    du = zeros(T, ndofs)
    ke = zeros(T, ndofs, ndofs)
    
    f!(fu, ke_temp, ue_interface, nms, ms)

    for i in 1:length(ue_interface)
        fill!(du, 0.0)
        fill!(f₊, 0.0)
        fill!(f₋, 0.0)
        du[i] = h
        f!(f₊, ke_temp, ue_interface + du, nms, deepcopy(ms))
        f!(f₋, ke_temp, ue_interface - du, nms, deepcopy(ms))
        ke[:,i] .= (f₊ - f₋)/(2h)
    end

    return fu, ke

end











