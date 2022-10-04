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


function Part{dim,T}(material::M, cellset::AbstractVector{Int}, element::E) where {dim,T,E,M}
    return Part{dim,T,E,M}(material, collect(cellset), element, PartCache{dim,T}(ndofs(element), Ferrite.nnodes(getcelltype(element))), PartVTKExport{dim,T}())
end

struct PartState{S<:AbstractElementState, M<:AbstractMaterialState} <: AbstractPartState
    elementstate::Vector{S}
    materialstates::Vector{M}
    stresses::Vector{SymmetricTensor{2,3,Float64,6}}
    strains::Vector{SymmetricTensor{2,3,Float64,6}}
end

get_fields(part::Part) = get_fields(part.element)
get_cellset(part::Part) = part.cellset

function construct_partstates(part::Part{dim,T,ET,MT}) where {dim,T,ET,MT}

    ncells = length(part.cellset)
    nqp = getnquadpoints(part.element)

    MaterialStateType = typeof( initial_material_state(part.material) )
    ElementStateType = elementstate_type(ET)

    states = Vector{PartState{ElementStateType,MaterialStateType}}(undef, ncells)

    for i in 1:ncells
        #@show MaterialStateType,ElementStateType
        _materialstates = Vector{MaterialStateType}(undef, nqp)
        _elementstates = Vector{ElementStateType}(undef, nqp)
        for j in 1:nqp
            _materialstates[j] = initial_material_state(part.material)
            _elementstates[j] = initial_element_state(part.element)
        end

        states[i] = PartState(_elementstates, _materialstates, zeros(SymmetricTensor{2,3,Float64,6}, nqp), zeros(SymmetricTensor{2,3,Float64,6}, nqp))
    end
    return states
end


function init_part!(part::Part, dh::Ferrite.AbstractDofHandler)
    celltype = typeof(dh.grid.cells[first(part.cellset)])
    vtk_celltype = Ferrite.cell_to_vtkcell(celltype)
    
    next_node_id = 1
    for cellid in part.cellset#CellIterator2(dh, part.element, part.cellset)
        cell = dh.grid.cells[cellid]

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

@enum ASSEMBLETYPE FORCEVEC STIFFMAT FSTAR DISSI

function assemble_stiffnessmatrix_and_forcevector!(dh::Ferrite.AbstractDofHandler, 
    part::Part,
    state::StateVariables)

    _assemble_part!(dh, part,state, STIFFMAT)

end

function assemble_forcevector!(dh::Ferrite.AbstractDofHandler, 
    part::Part,
    state::StateVariables)

    _assemble_part!(dh, part,state, FORCEVEC)

end

function assemble_fstar!(dh::Ferrite.AbstractDofHandler, 
    part::Part,
    state::StateVariables)

    _assemble_part!(dh, part,state, FSTAR)

end

function assemble_dissipation!(
    dh    ::Ferrite.AbstractDofHandler, 
    part  ::Part,
    state ::StateVariables)

    if !(is_dissipative(part.material) || is_dissipative(part.element))
        return 
    end

    _assemble_part!(dh, part, state, DISSI)

end

function _assemble_part!(dh::Ferrite.AbstractDofHandler, 
    part::Part{dim,T,ET,MT},
    state::StateVariables{T},
    assemtype::ASSEMBLETYPE) where {dim,T,ET,MT}

    assembler = start_assemble(state.system_arrays.Kⁱ, state.system_arrays.fⁱ, fillzero=false)
   
    ElementState = elementstate_type(ET)
    MaterialState = typeof( initial_material_state(part.material) )#materialstate_type(MT)

    #Extract variables
    element = part.element
    ue, Δue, due, ke, fe, coords, celldofs = (part.cache.ue, part.cache.Δue, part.cache.due, 
    part.cache.ke, part.cache.fe, part.cache.coords, part.cache.celldofs)

    Δt = state.Δt

    for (localid,cellid) in enumerate(part.cellset)
        
        partstate::PartState{ElementState, MaterialState} = state.partstates[cellid]

        materialstate = partstate.materialstates
        cellstate     = partstate.elementstate
        stresses      = partstate.stresses
        strains       = partstate.strains

        fill!(fe, 0.0)
        (assemtype == STIFFMAT) && fill!(ke, 0.0)

        Ferrite.cellcoords!(coords, dh, cellid)
        Ferrite.celldofs!(celldofs, dh, cellid)

        Δue .= state.v[celldofs] #Dont need both Δue and (v and Δt)
        ue .= state.d[celldofs]
        due .= state.v[celldofs]
        
        if assemtype == STIFFMAT
            integrate_forcevector_and_stiffnessmatrix!(element, cellstate, part.material, materialstate, stresses, strains, ke, fe, coords, Δue, ue, due, Δt)
            assemble!(assembler, celldofs, fe, ke)
        elseif assemtype == FORCEVEC
            integrate_forcevector!(element, cellstate, part.material, materialstate, fe, coords, Δue, ue, due, Δt)
            state.system_arrays.fⁱ[celldofs] += fe
        elseif assemtype == FSTAR
            error("Broken code, fix")
            prev_partstate::get_partstate_type(part) = state.prev_partstates[cellid]
            prev_materialstate = prev_partstate.materialstates

            integrate_fstar!(element, cellstate, part.material, prev_materialstate, fe, coords, Δue, ue, due, Δt)
            state.system_arrays.fᴬ[celldofs] += fe
        elseif assemtype == DISSI
            ge = Base.RefValue(zero(T))
            integrate_dissipation!(element, cellstate, part.material, materialstate, fe, ge, coords, Δue, ue, due, Δt)
            state.system_arrays.fᴬ[celldofs] += fe
            state.system_arrays.G[] += ge[]
        end

    end
    
end

function assemble_massmatrix!(dh::Ferrite.AbstractDofHandler, part::Part, state::StateVariables)

    assembler = start_assemble(state.system_arrays.M, fillzero=false)
    element = part.element

    dim = Ferrite.getdim(part)
    
    ue, due, me, fe, coords, celldofs = (part.cache.ue, part.cache.due, 
    part.cache.ke, part.cache.fe, part.cache.coords, part.cache.celldofs)

    for (localid,cellid) in enumerate(part.cellset)
        
        fill!(me, 0.0)

        partstate::get_partstate_type(part) = state.partstates[cellid]

        materialstate = partstate.materialstates
        cellstate     = partstate.elementstate

        Ferrite.cellcoords!(coords, dh, cellid)
        Ferrite.celldofs!(celldofs, dh, cellid)

        integrate_massmatrix!(element, cellstate, part.material, coords, me, ue, due)

        #assemble!(assembler, celldofs, me)
        state.system_arrays.M[celldofs, celldofs] += me
    end

end

function get_vtk_grid(dh::Ferrite.AbstractDofHandler, part::Part)
    return part.vtkexport.vtkcells, part.vtkexport.vtknodes
end

function post_part!(dh, part::Part, states::StateVariables)
    
end

function commit_part!(dh::Ferrite.AbstractDofHandler, part::Part, state::StateVariables)
    return nothing
end

#=
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
end=#

function get_vtk_field(dh::Ferrite.AbstractDofHandler, part::Part{dim,T}, state::StateVariables, field_name::Symbol) where {dim,T}
    fh = FieldHandler(get_fields(part), Set([1]))
    fieldidx = Ferrite.find_field(fh, field_name)
    @assert(fieldidx !== nothing)

    offset = Ferrite.field_offset(fh, field_name)
    fdim   = fh.fields[fieldidx].dim 

    n_vtk_nodes = length(part.vtkexport.vtknodes)
    #Special case displacement (it needs 3 datapoints)
    if field_name == :u
        data = zeros(T, (dim == 2 ? 3 : dim), n_vtk_nodes)
    else
        data = zeros(T, fdim, n_vtk_nodes)
    end
    _get_vtk_field!(data, dh, part, state, offset, fdim)
    return data
end

function _get_vtk_field!(data::Matrix, dh::Ferrite.AbstractDofHandler, part::Part{dim,T}, state::StateVariables, offset::Int, nvars::Int) where {dim,T}

    celldofs = part.cache.celldofs
    for cellid in part.cellset
        cell = dh.grid.cells[cellid]
        
        celldofs!(celldofs, dh, cellid)
        ue = state.d[celldofs]
        counter = 1
        for (i,nodeid) in enumerate(cell.nodes)
            local_id = part.vtkexport.nodeid_mapper[nodeid]
            for d in 1:nvars
                data[d, local_id] = ue[counter + offset]
                counter += 1
            end
        end
    end
    
end


function get_vtk_celldata(part::Part{dim}, output::VTKCellOutput{<:StressOutput}, state::StateVariables{T}, globaldata) where {dim,T}
    
    npartcells = length(part.cellset)
    M = 9 #Always symmetric 3d stress tensor
    data = Matrix{T}(undef, M, npartcells)
    
    for (ic, cellid) in enumerate(part.cellset)
        stresses = state.partstates[cellid].stresses
        data[:, ic] .= stresses |> output.func |> (x) -> reinterpret(T, x) |> collect |> vec
    end

    return data
end

function get_vtk_celldata(part::Part, output::VTKCellOutput{<:MaterialStateOutput}, state::StateVariables{T}, globaldata) where T

    _cellid = first(part.cellset)
    first_state = first(state.partstates[_cellid].materialstates)
    if !hasproperty(first_state, output.type.field) 
        return nothing
    end
    
    StateType = typeof(first_state)
    npartcells = length(part.cellset)
    ncomp = length(getproperty(first_state, output.type.field))
    data = Matrix{T}(undef, ncomp, npartcells)
    
    for (ic, cellid) in enumerate(part.cellset)
        matstates = state.partstates[cellid].materialstates
        data[:, ic] .= getproperty.(matstates, output.type.field) |> output.func |> (x) -> reinterpret(T, x) |> collect |> vec
    end

    return data
end

function get_vtk_nodedata(part::Part{dim}, output::VTKNodeOutput{<:StressOutput}, state::StateVariables{T}, globaldata) where {dim,T}
   
    #Extract stresses to interpolate
    data = Vector{SymmetricTensor{2,3,T,6}}[]
    for (ic, cellid) in enumerate(part.cellset)
        stresses = state.partstates[cellid].stresses
        push!(data, stresses)
    end

    _get_vtk_nodedata(part, data, globaldata)
end

function get_vtk_nodedata(part::Part{dim}, output::VTKNodeOutput{<:MaterialStateOutput}, state::StateVariables{T}, globaldata) where {dim,T}

    #TODO: Move this check before simulatoin?
    _cellid = first(part.cellset)
    first_state = first(state.partstates[_cellid].materialstates)
    if !hasproperty(first_state, output.type.field) 
        return nothing
    end
    first_field = getproperty(first_state, output.type.field)
    FieldDataType = typeof(first_field)
  
    #Extract field to interpolate
    data = Vector{FieldDataType}[]
    for (ic, cellid) in enumerate(part.cellset)
        matstates = state.partstates[cellid].materialstates
        field_states = getproperty.(matstates, output.type.field)
        push!(data, field_states)
    end

    _get_vtk_nodedata(part, data, globaldata)
end

function _get_vtk_nodedata(part::Part{dim}, data::Vector{Vector{FieldDataType}}, globaldata) where {dim,FieldDataType}
    
    #Set up quadrature rule
    celltype = getcelltype(part.element)
    geom_ip = Ferrite.default_interpolation(celltype)
    qr = getquadraturerule(part.element)
    #refshape = Ferrite.getrefshape(geom_ip)
    #nqp = length(first(data))
    #qp_order = convert(Int, nqp^(1/dim))
    #qr = QuadratureRule{dim, refshape}(qp_order) #Does not work for shells?

    #Project
    #= Projection for cohseive zone elements
        L2Projector assumes CellScalarValues, But Should the depracted methods works with coustum SurfaceVectorValues...
    if typeof(geom_ip) <: CohesiveZoneInterpolation
        fe_values = part.element.cv
        projector = Ferrite.L2Projector(fe_values, geom_ip, globaldata.grid, part.cellset)
        n = getnbasefunctions(fe_values)
        @show n
        #@show ndofs_per_cell(globaldata.dh, first(part.cellset))
        data_nodes = project(projector, data; project_to_nodes=true); # TODO: this should be default.
    end=#

    projector = Ferrite.L2Projector(geom_ip, globaldata.grid; set = part.cellset)
    data_nodes = project(projector, data, qr; project_to_nodes=true); # TODO: this should be default.

    #Reorder to the parts vtk
    nvtknodes = length(part.vtkexport.nodeid_mapper)
    vtk_node_data = zeros(FieldDataType, nvtknodes)# Matrix{T}(undef, ncomp, nvtknodes)
    for (ic,cellid) in enumerate(part.cellset)
        for nodeid in globaldata.grid.cells[cellid].nodes
            vtknodeid = part.vtkexport.nodeid_mapper[nodeid]
            vtk_node_data[vtknodeid] = data_nodes[nodeid]
        end
    end

    return vtk_node_data
end

function get_vtk_celldata(dh::Ferrite.AbstractDofHandler, part::Part, state::StateVariables) where {dim_p,dim_s,T}
    return nothing, nothing
end

function get_vtk_nodedata(dh::Ferrite.AbstractDofHandler, part::Part, state::StateVariables) where {dim_p,dim_s,T}
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

function build_output!(part::Part, output::OutputData, state::StateVariables, globaldata)

    

end