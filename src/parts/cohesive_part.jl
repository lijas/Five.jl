
#Cohesive part, just override some plotting stuff
function init_part!(part::Part{dim,T,<:CohesiveElement}, dh::Ferrite.AbstractDofHandler) where {dim,T}
    @assert(dim==2)
    celltype = typeof(dh.grid.cells[first(part.cellset)])
    vtk_celltype = Ferrite.VTKCellTypes.VTK_LINE
    @assert(vtk_celltype.nodes == 2)

    next_node_id = 1
    for cellid in part.cellset#CellIterator2(dh, part.element, part.cellset)
        cell = dh.grid.cells[cellid]
        new_ids = Int[]
        for nodeid in cell.nodes[1:2]
            if !haskey(part.vtkexport.nodeid_mapper, nodeid)
                part.vtkexport.nodeid_mapper[nodeid] = next_node_id
                push!(new_ids, next_node_id)
                next_node_id += 1
                push!(part.vtkexport.vtknodes, dh.grid.nodes[nodeid].x)
            else
                _new_id = part.vtkexport.nodeid_mapper[nodeid]
                push!(new_ids, _new_id)
            end
        end
        push!(part.vtkexport.vtkcells, MeshCell(vtk_celltype, new_ids))
    end

    resize!(part.cache.coords, Ferrite.nnodes(celltype))
end

function get_vtk_displacements(dh::Ferrite.AbstractDofHandler, part::Part{dim,T,<:CohesiveElement}, state::StateVariables) where {dim,T}
    @assert(length(get_fields(part.element)) == 1 && get_fields(part.element)[1].name == :u)

    node_coords = zeros(Vec{dim,T}, length(part.vtkexport.vtknodes))

    celldofs = part.cache.celldofs

    for cellid in part.cellset
        cell = dh.grid.cells[cellid]

        celldofs!(celldofs, dh, cellid)
        ue = state.d[celldofs]
        uvec = reinterpret(Vec{dim,T}, ue[1:4*dim])

        nodeid = cell.nodes[1]
        local_id = part.vtkexport.nodeid_mapper[nodeid]
        node_coords[local_id] = (uvec[1] + uvec[3])/2

        nodeid = cell.nodes[2]
        local_id = part.vtkexport.nodeid_mapper[nodeid]
        node_coords[local_id] = (uvec[2] + uvec[4])/2
    end
    
    return node_coords
end