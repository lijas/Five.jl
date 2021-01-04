export Part
export IGAPart, PartState

struct PartCache{dim,T}
    ue::Vector{T}
    due::Vector{T}
    Δue::Vector{T}
    fe::Vector{T}
    ke::Matrix{T}
    celldofs::Vector{Int}
    coords::Vector{Vec{dim,T}}
end

struct PartVTKExport{dim,T}
    nodeid_mapper::Dict{Int,Int}
    vtkcells::Vector{MeshCell}
    vtknodes::Vector{Vec{dim,T}}
end

PartVTKExport{dim,T}() where {dim,T} = PartVTKExport{dim,T}(Dict{Int,Int}(), MeshCell[], Vec{dim,T}[])
PartCache{dim,T}(ndofs, ncoords) where {dim,T} = PartCache{dim,T}(zeros(T,ndofs), zeros(T,ndofs), zeros(T,ndofs), zeros(T,ndofs), zeros(T,ndofs,ndofs), zeros(Int,ndofs), zeros(Vec{dim,T},ncoords))

#Name it Part instead of Fe-part because it is the standard...
struct Part{dim, T, E<:AbstractElement, M<:AbstractMaterial} <: AbstractPart{dim}
    #id::Int
    material::M
    cellset::Vector{Int}
    element::E

    #For vtk-ploting
    cache::PartCache{dim,T}
    vtkexport::PartVTKExport{dim,T}
end

function Part{dim,T}(; 
    material::AbstractMaterial,
    element::AbstractElement,
    cellset::Vector{Int}
    ) where {dim,T}

   return Part{dim,T}(material, cellset, element)
end

struct IGAPart{dim, T, E<:AbstractElement, M<:AbstractMaterial} <: AbstractPart{dim}
 
    material::M
    cellset::Vector{Int}
    element::E

    Cb::Vector{IGA.BezierExtractionOperator{T}}
    cv_plot::CellScalarValues 
end

const FEPart{dim} = Union{IGAPart{dim}, Part{dim}}# where dim

function Part{dim,T}(material::M, cellset::AbstractVector{Int}, element::E) where {dim,T,E,M}
    return Part{dim,T,E,M}(material, collect(cellset), element, PartCache{dim,T}(ndofs(element), getncoords(element)), PartVTKExport{dim,T}())
end

struct PartState{S<:AbstractElementState, M<:AbstractMaterialState} <: AbstractPartState
    elementstate::S 
    materialstates::Vector{M} 
end

getmaterialstate(partstate::PartState, field::Symbol) =  getproperty.(partstate.materialstates, field)

get_partstate_type(part::Part{dim,T,E,M}) where {dim,T,E,M} = return PartState{get_elementstate_type(part.element), get_material_state_type(part.material)}
get_partstate_type(part::IGAPart{dim,T,E,M}) where {dim,T,E,M} = return PartState{get_elementstate_type(part.element), get_material_state_type(part.material)}

get_fields(part::FEPart) = get_fields(part.element)
get_cellset(part::FEPart) = part.cellset

function construct_partstates(part::FEPart)

    ncells = length(part.cellset)
    nqp = getnquadpoints(part.element)

    materialtype = get_material_state_type(part.material)
    ElementStateType = get_elementstate_type(part.element)

    states = Vector{PartState{ElementStateType,materialtype}}(undef, ncells)

    for i in 1:ncells
        _cellstate = ElementStateType(part.element)

        _materialstates = Vector{materialtype}(undef, nqp)
        for j in 1:nqp
            _materialstates[j] = materialtype(part.material)
        end

        states[i] = PartState(_cellstate, _materialstates)
    end
    return states
end


function init!(part::Part, dh::JuAFEM.AbstractDofHandler)
    celltype = JuAFEM.cell_to_vtkcell(getcelltype(part.element))
    
    #cls = MeshCell[]
    #nodeid_mapper = Dict{Int,Int}()
    #node_coords = Vec{dim,T}[]
    next_node_id = 1
    for cell in CellIterator2(dh, part.element, part.cellset)
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
        push!(part.vtkexport.vtkcells, MeshCell(celltype, copy(new_ids[1:celltype.nodes])))
    end
     #@show typeof(part.element)
     #@show part.vtknode_coords
end

@enum ASSEMBLETYPE FORCEVEC STIFFMAT FSTAR

function assemble_stiffnessmatrix_and_forcevector!(dh::JuAFEM.AbstractDofHandler, 
    part::FEPart,
    state::StateVariables)

    _assemble_part!(dh, part,state, STIFFMAT)

end

function assemble_forcevector!(dh::JuAFEM.AbstractDofHandler, 
    part::FEPart,
    state::StateVariables)

    _assemble_part!(dh, part,state, FORCEVEC)

end

function assemble_fstar!(dh::JuAFEM.AbstractDofHandler, 
    part::FEPart,
    state::StateVariables)

    _assemble_part!(dh, part,state, FSTAR)

end

function _assemble_part!(dh::JuAFEM.AbstractDofHandler, 
    part::FEPart{dim},
    state::StateVariables{T},
    assemtype::ASSEMBLETYPE) where {dim,T}

    assembler = start_assemble(state.system_arrays.Kⁱ, state.system_arrays.fⁱ, fillzero=false)
    
    #Extract variables
    element = part.element
    ue, Δue, due, ke, fe, coords, celldofs = (part.cache.ue, part.cache.Δue, part.cache.due, 
    part.cache.ke, part.cache.fe, part.cache.coords, part.cache.celldofs)

    Δt = state.Δt

    for (localid,cellid) in enumerate(part.cellset)
        
        partstate::get_partstate_type(part) = state.partstates[cellid]
        ⁿpartstate::get_partstate_type(part) = state.prev_partstates[cellid]

        materialstate      = partstate.materialstates
        cellstate          = partstate.elementstate

        fill!(fe, 0.0)
        (assemtype == STIFFMAT) && fill!(ke, 0.0)

        JuAFEM.cellcoords!(coords, dh, cellid)
        JuAFEM.celldofs!(celldofs, dh, cellid)

        Δue .= state.Δd[celldofs]
        ue .= state.d[celldofs]
        due .= state.v[celldofs]

        if typeof(part) <: IGAPart
            Ce = get_bezier_operator(part, localid)
            IGA.set_bezier_operator!(part.element.cv, Ce)
            coords .= IGA.compute_bezier_points(Ce, coords)
        end
        
        if assemtype == STIFFMAT
            integrate_forcevector_and_stiffnessmatrix!(element, cellstate, part.material, materialstate, ke, fe, coords, Δue, ue, due, Δt)
            assemble!(assembler, celldofs, fe, ke)
        elseif assemtype == FORCEVEC
            integrate_forcevector!(element, cellstate, part.material, materialstate, fe, coords, Δue, ue, due, Δt)
            state.system_arrays.fⁱ[celldofs] += fe
        elseif assemtype == FSTAR
            integrate_fstar!(element, cellstate, part.material, materialstate, fe, coords, Δue, ue, due, Δt)
            state.system_arrays.fⁱ[celldofs] += fe
        end

    end
    
end

function  assemble_massmatrix!(dh::JuAFEM.AbstractDofHandler, part::FEPart, state::StateVariables) where T

    assembler = start_assemble(state.system_arrays.M, fillzero=false)
    element = part.element

    dim = JuAFEM.getdim(part)
    
    #preallocate stuff (could be stored in Part)
    me = zeros(T, ndofs(element), ndofs(element))
    ue = zeros(T, ndofs(element))
    due = zeros(T, ndofs(element))

    coords = zeros(Vec{dim,T}, JuAFEM.nnodes_per_cell(dh, first(part.cellset)))
    celldofs = zeros(Int, ndofs(element))

    for (localid,celldata) in enumerate(CellIterator(dh,part.cellset))
        
        fill!(me, 0.0)

        JuAFEM.cellcoords!(coords, dh, cellid(celldata))
        JuAFEM.celldofs!(celldofs, dh, cellid(celldata))

        if typeof(part) <: IGAPart
            Ce = get_bezier_operator(part, localid)
            IGA.set_bezier_operator!(part.element.cv, Ce)
            coords .= IGA.compute_bezier_points(Ce, coords)
        end

        integrate_massmatrix!(element, get_elementstate_type(element)(), part.material, coords, me, ue, due)

        assemble!(assembler, celldofs, me)
    end

end

function get_vtk_grid(dh::JuAFEM.AbstractDofHandler, part::Part)
    return part.vtkexport.vtkcells, part.vtkexport.vtknodes
end

function post_part!(dh, part::FEPart, states::StateVariables)
    
end

function commit_part!(dh::JuAFEM.AbstractDofHandler, part::FEPart, state::StateVariables)
    return nothing
end

function get_vtk_displacements(dh::JuAFEM.AbstractDofHandler, part::Part{dim,T}, state::StateVariables) where {dim,T}
    @assert(length(get_fields(part.element)) == 1 && get_fields(part.element)[1].name == :u)

    node_coords = zeros(Vec{dim,T}, length(part.vtkexport.vtknodes))

    for cell in CellIterator2(dh, part.element, part.cellset)
        ue = state.d[cell.celldofs]
        ue_vec = reinterpret(Vec{dim,T}, ue)

        coords = ue_vec

        for (i,nodeid) in enumerate(cell.nodes)
            local_id = part.vtkexport.nodeid_mapper[nodeid]
            node_coords[local_id] = coords[i]
        end
    end
    
    return node_coords
end

function init!(part::Part{dim,T,<:CohesiveElement}, dh::MixedDofHandler{dim,C,T}) where {dim,C,T}
    
    @assert(dim==2)
    celltype = JuAFEM.VTKCellTypes.VTK_LINE
    
    next_node_id = 1
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
        
    end
end

function get_vtk_part_grid(part::Part{dim,T,<:CohesiveElement}, dh::MixedDofHandler{dim,T}) where {dim,T}    
    return part.vtkexport.vtkcells, part.vtkexport.vtknodes
end

function get_vtk_part_point_data(part::Part{dim,T,<:CohesiveElement}, dh::MixedDofHandler{dim,T}, di) where {dim,T}
    @assert(length(part.element.fields) == 1 && part.element.fields[1].name == :u)

    node_coords = zeros(Vec{dim,T}, length(part.vtkexport.vtknodes))

    for cell in CellIterator(dh, part.element, part.cellset)
        ue = di[cell.celldofs]
        ue_vec = reinterpret(Vec{dim,T}, ue)

        coords = ue_vec

        for (i,nodeid) in enumerate(cell.nodes)
            local_id = part.vtkexport.nodeid_mapper[nodeid]
            node_coords[local_id] = coords[i]
        end
    end
    
    return node_coords
end

function get_vtk_celldata(part::FEPart, output::VTKCellOutput{<:MaterialStateOutput}, state::StateVariables{T}, globaldata) where T

    first_state = first(first(state.partstates).materialstates)
    if !hasproperty(first_state, output.type.field) 
        return nothing
    end
    
    StateType = typeof(first_state)
    npartcells = length(part.cellset)
    ncomp = length(getproperty(first_state, output.type.field))
    data = Matrix{T}(undef, npartcells, ncomp)
    
    for (ic, cellid) in enumerate(part.cellset)
        matstates = state.partstates[cellid].materialstates
        data[ic, :] = getproperty.(matstates, output.type.field) |> output.func |> (x) -> reinterpret(T, x) |> vec
    end

    return data
end

function get_vtk_nodedata(part::FEPart, output::VTKNodeOutput{<:MaterialStateOutput}, state::StateVariables{T}, globaldata) where T

    

    return data
end

function get_vtk_celldata(dh::JuAFEM.AbstractDofHandler, part::FEPart, state::StateVariables) where {dim_p,dim_s,T}
    return nothing, nothing
end

function get_vtk_nodedata(dh::JuAFEM.AbstractDofHandler, part::FEPart, state::StateVariables) where {dim_p,dim_s,T}
    return nothing, nothing
end


#=function init!(part::Part{<:SolidElement}, dh::DofHandler{dim,T}) where {dim,T}
    
    @assert(dim==2)
end

function get_vtk_part_grid(part::Part{<:SolidElement}, dh::DofHandler{dim,T}) where {dim,T}    
    return part.vtkcells, part.vtknode_coords
end

function get_vtk_part_point_data(part::Part{<:SolidElement}, dh::DofHandler{dim,T}, di) where {dim,T}
    return []
end=#