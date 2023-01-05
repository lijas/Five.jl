
#Cohesive part, just override some plotting stuff
function init_part!(part::Part{dim,T,<:CohesiveElement}, dh::Ferrite.AbstractDofHandler) where {dim,T}
   #= celltype = typeof(dh.grid.cells[first(part.cellset)])
    vtk_celltype = Ferrite.cell_to_vtkcell(celltype)

    next_node_id = 1
    for cellid in part.cellset#CellIterator2(dh, part.element, part.cellset)
        cell = dh.grid.cells[cellid]
        new_ids = Int[]
        for i in 1:(Ferrite.nnodes(celltype) ÷ 2) # Mid surface
            nodeid = cell.nodes[i]
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
 #       @show part.vtkexport.vtknodes[new_ids]
        push!(part.vtkexport.vtkcells, MeshCell(vtk_celltype, new_ids))
    end
#asdf
    resize!(part.cache.coords, Ferrite.nnodes(celltype))=#
end

function get_vtk_celldata(part::Part{dim,T,<:CohesiveElement}, output::VTKCellOutput{<:StressOutput}, state::StateVariables{T}, globaldata) where {dim,T}
    return nothing
end

function get_vtk_nodedata(part::Part{dim,T,<:CohesiveElement}, output::VTKNodeOutput{<:StressOutput}, state::StateVariables{T}, globaldata) where {dim,T}
    return nothing
end

function get_vtk_displacements(dh::Ferrite.AbstractDofHandler, part::Part{dim,T,<:CohesiveElement}, state::StateVariables) where {dim,T}
    @assert(length(get_fields(part.element)) == 1 && get_fields(part.element)[1].name == :u)

    node_disp = zeros(Vec{dim,T}, length(part.vtkexport.vtknodes))

    cell = dh.grid.cells[first(part.cellset)]
    celltype = typeof(cell)
    n_mid_surface_nodes = (Ferrite.nnodes(celltype) ÷ 2)

    celldofs = part.cache.celldofs

    for cellid in part.cellset
        cell = dh.grid.cells[cellid]

        celldofs!(celldofs, dh, cellid)
        ue = state.d[celldofs]
        uvec = reinterpret(Vec{dim,T}, ue)

        for i in 1:n_mid_surface_nodes
            nodeid = cell.nodes[i]
            local_id = part.vtkexport.nodeid_mapper[nodeid] 
            node_disp[local_id] = (uvec[i] + uvec[i+n_mid_surface_nodes])/2 
        end
    end
    
    return node_disp
end


function _get_vtk_field!(data::Matrix, dh::Ferrite.AbstractDofHandler, part::Part{dim,T,<:CohesiveElement}, state::StateVariables, offset::Int, nvars::Int) where {dim,T}

    celldofs = part.cache.celldofs

    for cellid in part.cellset
        cell = dh.grid.cells[cellid]
        
        celldofs!(celldofs, dh, cellid)
        ue = state.d[celldofs]
        counter = 1
        for i in 1:(Ferrite.nnodes(cell) ÷ 2) # Mid surface
            nodeid = cell.nodes[i]
            local_id = part.vtkexport.nodeid_mapper[nodeid]
            for d in 1:nvars
                data[d, local_id] = ue[counter + offset]
                counter += 1
            end
        end
    end
    
end

function collect_nodedata!(data::Vector{FT}, part::Part{dim,T,<:CohesiveElement}, output::StressOutput, state::StateVariables{T}, globaldata) where {dim,FT,T} 
end