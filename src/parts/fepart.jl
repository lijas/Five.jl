export Part
export PartState

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
    cellset::AbstractVector{Int}
    ) where {dim,T}

   return Part{dim,T}(material, collect(cellset), element)
end

#Not implemented:
struct IGAPart{dim, T, E<:AbstractElement, M<:AbstractMaterial} <: AbstractPart{dim}
 
    material::M
    cellset::Vector{Int}
    element::E

    #Cb::Vector{IGA.BezierExtractionOperator{T}}
    #cv_plot::CellScalarValues 
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
            _materialstates[j] = getmaterialstate(part.material)
        end

        states[i] = PartState(_cellstate, _materialstates)
    end
    return states
end


function init_part!(part::Part, dh::Ferrite.AbstractDofHandler)
    celltype = typeof(dh.grid.cells[first(part.cellset)])
    
    next_node_id = 1
    for cellid in part.cellset#CellIterator2(dh, part.element, part.cellset)
        cell = dh.grid.cells[cellid]
        vtk_celltype = Ferrite.cell_to_vtkcell(celltype)

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
        push!(part.vtkexport.vtkcells, MeshCell(vtk_celltype, copy(new_ids[1:vtk_celltype.nodes])))
    end


    resize!(part.cache.coords, Ferrite.nnodes(celltype))
end

@enum ASSEMBLETYPE FORCEVEC STIFFMAT FSTAR DISSI ELASTIC

function assemble_stiffnessmatrix_and_forcevector!(dh::Ferrite.AbstractDofHandler, 
    part::FEPart,
    state::StateVariables)

    _assemble_part!(dh, part,state, STIFFMAT)

end

function assemble_forcevector!(dh::Ferrite.AbstractDofHandler, 
    part::FEPart,
    state::StateVariables)

    _assemble_part!(dh, part,state, FORCEVEC)

end

function assemble_fstar!(dh::Ferrite.AbstractDofHandler, 
    part::FEPart,
    state::StateVariables)

    _assemble_part!(dh, part,state, FSTAR)

end

function assemble_dissipation!(dh::Ferrite.AbstractDofHandler, 
    part::FEPart,
    state::StateVariables)

    if !(part.material |> is_dissipative)
        return 
    end

    _assemble_part!(dh, part, state, DISSI)

end

function assemble_elastic!(dh::Ferrite.AbstractDofHandler, 
    part::FEPart,
    state::StateVariables)

    _assemble_part!(dh, part, state, ELASTIC)

end

function _assemble_part!(dh::Ferrite.AbstractDofHandler, 
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

        materialstate = partstate.materialstates
        cellstate     = partstate.elementstate

        fill!(fe, 0.0)
        (assemtype == STIFFMAT) && fill!(ke, 0.0)
        #@show typeof(element)
        #@show length(celldofs)
        #@show ndofs_per_cell(dh,cellid)
        Ferrite.cellcoords!(coords, dh, cellid)
        Ferrite.celldofs!(celldofs, dh, cellid)

        Δue .= state.Δd[celldofs]
        ue .= state.d[celldofs]
        due .= state.v[celldofs]
        
        if assemtype == STIFFMAT
            integrate_forcevector_and_stiffnessmatrix!(element, cellstate, part.material, materialstate, ke, fe, coords, Δue, ue, due, Δt)
            assemble!(assembler, celldofs, fe, ke)
        elseif assemtype == FORCEVEC
            integrate_forcevector!(element, cellstate, part.material, materialstate, fe, coords, Δue, ue, due, Δt)
            state.system_arrays.fⁱ[celldofs] += fe
        elseif assemtype == FSTAR
            prev_partstate::get_partstate_type(part) = state.prev_partstates[cellid]
            prev_materialstate = prev_partstate.materialstates

            integrate_fstar!(element, cellstate, part.material, prev_materialstate, fe, coords, Δue, ue, due, Δt)
            state.system_arrays.fᴬ[celldofs] += fe
        elseif assemtype == DISSI
            prev_partstate = state.prev_partstates[cellid]
            prev_materialstate = prev_partstate.materialstates

            ge = Base.RefValue(zero(T))
            integrate_dissipation!(element, cellstate, part.material, prev_materialstate, fe, ge, coords, Δue, ue, due, Δt)
            state.system_arrays.fᴬ[celldofs] += fe
            state.system_arrays.G[] += ge[]
        elseif assemtype == ELASTIC
            prev_partstate = state.prev_partstates[cellid]
            prev_materialstate = prev_partstate.materialstates

            ge = Base.RefValue(zero(T))
            integrate_elastic!(element, cellstate, part.material, prev_materialstate, fe, ge, coords, Δue, ue, due, Δt)
            state.system_arrays.fᴬ[celldofs] += fe
            state.system_arrays.G[] += ge[]
        end

    end
    
end

function assemble_massmatrix!(dh::Ferrite.AbstractDofHandler, part::FEPart, state::StateVariables) where T

    assembler = start_assemble(state.system_arrays.M, fillzero=false)
    element = part.element

    dim = Ferrite.getdim(part)
    
    #preallocate stuff (could be stored in Part)
    me = zeros(T, ndofs(element), ndofs(element))
    ue = zeros(T, ndofs(element))
    due = zeros(T, ndofs(element))

    coords = zeros(Vec{dim,T}, Ferrite.nnodes_per_cell(dh, first(part.cellset)))
    celldofs = zeros(Int, ndofs(element))

    for (localid,celldata) in enumerate(CellIterator(dh,part.cellset))
        
        fill!(me, 0.0)

        Ferrite.cellcoords!(coords, dh, cellid(celldata))
        Ferrite.celldofs!(celldofs, dh, cellid(celldata))

        integrate_massmatrix!(element, get_elementstate_type(element)(), part.material, coords, me, ue, due)

        assemble!(assembler, celldofs, me)
    end

end

function get_vtk_grid(dh::Ferrite.AbstractDofHandler, part::Part)
    return part.vtkexport.vtkcells, part.vtkexport.vtknodes
end

function post_part!(dh, part::FEPart, states::StateVariables)
    
end

function commit_part!(dh::Ferrite.AbstractDofHandler, part::FEPart, state::StateVariables)
    return nothing
end

function get_vtk_displacements(dh::Ferrite.AbstractDofHandler, part::Part{dim,T}, state::StateVariables) where {dim,T}
    @assert(length(get_fields(part.element)) == 1 && get_fields(part.element)[1].name == :u)

    node_coords = zeros(Vec{dim,T}, length(part.vtkexport.vtknodes))

    celldofs = part.cache.celldofs

    for cellid in part.cellset#CellIterator2(dh, part.element, part.cellset)
        cell = dh.grid.cells[cellid]

        celldofs!(celldofs, dh, cellid)
        ue = state.d[celldofs]
        ue_vec = reinterpret(Vec{dim,T}, ue)
        for (i,nodeid) in enumerate(cell.nodes)
            local_id = part.vtkexport.nodeid_mapper[nodeid]
            node_coords[local_id] = ue_vec[i]
        end
    end
    
    return node_coords
end

function get_vtk_celldata(part::FEPart, output::VTKCellOutput{<:MaterialStateOutput}, state::StateVariables{T}, globaldata) where T

    _cellid = first(part.cellset)
    first_state = first(state.partstates[_cellid].materialstates)
    if !hasproperty(first_state, output.type.field) 
        return nothing
    end
    
    StateType = typeof(first_state)
    npartcells = length(part.cellset)
    ncomp = length(getproperty(first_state, output.type.field))
    data = Matrix{T}(undef, npartcells, ncomp)
    
    for (ic, cellid) in enumerate(part.cellset)
        matstates = state.partstates[cellid].materialstates
        data[ic, :] .= getproperty.(matstates, output.type.field) |> output.func |> (x) -> reinterpret(T, x) |> collect |> vec
    end

    return data
end

function get_vtk_nodedata(part::FEPart{dim}, output::VTKNodeOutput{<:MaterialStateOutput}, state::StateVariables{T}, globaldata) where {dim,T}
    
    _cellid = first(part.cellset)
    first_state = first(state.partstates[_cellid].materialstates)
    if !hasproperty(first_state, output.type.field) 
        return nothing
    end
    first_field = getproperty(first_state, output.type.field)
    
    ncomp = length(first_field)
    FieldDataType = typeof(first_field)
    #celltype = getcelltype(part.element)
    celltype = typeof(globaldata.grid.cells[first(part.cellset)])
    geom_ip = Ferrite.default_interpolation(celltype)
    refshape = Ferrite.getrefshape(geom_ip)# RefCube#getrefshape(geom_ip)
    nqp = length(state.partstates[_cellid].materialstates)
    #qp_order = convert(Int, nqp^(1/dim))
    qr = QuadratureRule{dim, refshape}(2)#qp_order)
    cellvalues = CellScalarValues(qr, geom_ip)
    projector = L2Projector(cellvalues, geom_ip, globaldata.grid, part.cellset)

    #Extract field to interpolate
    data = Vector{FieldDataType==Float64 ? Vec{1,T} : FieldDataType}[]
    for (ic, cellid) in enumerate(part.cellset)
        matstates = state.partstates[cellid].materialstates
        field_states = getproperty.(matstates, output.type.field)
        if FieldDataType == Float64
            push!(data, reinterpret(Vec{1,T}, field_states))
        else
            push!(data, field_states)
        end
    end

    data_nodes = project(data, projector)
    nvtknodes = length(part.vtkexport.nodeid_mapper)
    vtk_node_data = Matrix{T}(undef, ncomp, nvtknodes)
    for (ic,cellid) in enumerate(part.cellset)
        for nodeid in globaldata.grid.cells[cellid].nodes
            vtknodeid = part.vtkexport.nodeid_mapper[nodeid]
            vtk_node_data[:,vtknodeid] .= reinterpret(T, data_nodes[nodeid])
        end
    end
    return vtk_node_data
end

function get_vtk_celldata(dh::Ferrite.AbstractDofHandler, part::FEPart, state::StateVariables) where {dim_p,dim_s,T}
    return nothing, nothing
end

function get_vtk_nodedata(dh::Ferrite.AbstractDofHandler, part::FEPart, state::StateVariables) where {dim_p,dim_s,T}
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