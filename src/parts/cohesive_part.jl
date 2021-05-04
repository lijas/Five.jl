
#Cohesive part, just override some plotting stuff

function init_part!(part::Part{dim,T,<:CohesiveElement}, dh::Ferrite.AbstractDofHandler) where {dim,T}
    return nothing
    #=@assert(dim==2)
    celltype = Ferrite.VTKCellTypes.VTK_LINE
    
    next_node_id = 1
    @show part.cellset
    for cell in CellIterator(dh, part.cellset)
        new_ids = Int[]
        for nodeid in cell.nodes
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
        #@show cell.nodes
        push!(part.vtkexport.vtkcells, MeshCell(celltype, new_ids[1:2]))
        push!(part.vtkexport.vtkcells, MeshCell(celltype, new_ids[2:3]))
        push!(part.vtkexport.vtkcells, MeshCell(celltype, new_ids[3:4]))
        push!(part.vtkexport.vtkcells, MeshCell(celltype, new_ids[[4,1]]))
        
    end=#
end

function get_vtk_displacements(dh::Ferrite.AbstractDofHandler, part::Part{dim,T,<:CohesiveElement}, state::StateVariables) where {dim,T}
    return Vec{dim,T}[]
end

