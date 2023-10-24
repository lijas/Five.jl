
struct PartCache{dim,T,E}
    ue::Vector{T}
    due::Vector{T}
    Δue::Vector{T}
    fe::Vector{T}
    ke::Matrix{T}
    celldofs::Vector{Int}
    coords::Vector{Vec{dim,T}}
    element::E
end

function PartCache{dim,T}(ndofs, ncoords, element::E) where {dim,T,E} 
    PartCache{dim,T,E}(
        zeros(T,ndofs), 
        zeros(T,ndofs), 
        zeros(T,ndofs), 
        zeros(T,ndofs), 
        zeros(T,ndofs,ndofs), 
        zeros(Int,ndofs), 
        zeros(Vec{dim,T},ncoords),
        deepcopy(element)
    )
end

#Name it Part instead of Fe-part because it is the standard...
struct Part{dim, T, E<:AbstractElement, M<:AbstractMaterial} <: AbstractPart{dim}
    material::M
    cellset::Vector{Int}
    threadsets::Vector{Vector{Int}}
    element::E

    cache::Vector{PartCache{dim,T,E}} #One for each thread
    geometry::Optional{Any} # TODO: Where should the visualization geometry live?
end

function Part(; 
    material::M,
    element::E,
    cellset,
    geometry=nothing
    ) where {E,M}

    dim = Ferrite.getdim(element)
    T = Float64
    _set = collect(cellset)
    sort!(_set) # YOLO
    return Part{dim,T,E,M}(
        material, 
        _set,
        Vector{Int}[],
        element, 
        PartCache{dim,T}[],
        geometry)
end

struct PartState{ES<:AbstractElementState, MS<:AbstractMaterialState} <: AbstractPartState
    elementstate   ::Vector{Vector{ES}}
    materialstates ::Vector{Vector{MS}}
    stresses       ::Vector{Vector{SymmetricTensor{2,3,Float64,6}}}
    strains        ::Vector{Vector{SymmetricTensor{2,3,Float64,6}}}
end

get_fields(part::Part) = get_fields(part.element)
get_cellset(part::Part) = part.cellset

function construct_partstates(part::Part{dim,T,ET,MT}) where {dim,T,ET,MT}

    ncells = length(part.cellset)
    nqp    = getnquadpoints(part.element)

    materialstate_type = get_material_state_type(MT)
    elementstate_type = get_element_state_type(ET)

    estates = Vector{elementstate_type}[]
    mstates = Vector{materialstate_type}[]
    stresses = Vector{SymmetricTensor{2,3,Float64,6}}[]
    strains = Vector{SymmetricTensor{2,3,Float64,6}}[]
    for i in 1:ncells
        push!(estates, [initial_element_state(part.element) for i in 1:nqp])
        push!(mstates, [initial_material_state(part.material) for i in 1:nqp])
        push!(stresses, zeros(SymmetricTensor{2,3,Float64,6}, nqp))
        push!(strains, zeros(SymmetricTensor{2,3,Float64,6}, nqp))
    end

    pstate = PartState(estates, mstates, stresses, strains)
    return pstate
end

function init_part!(part::Part{dim, T}, dh::Ferrite.AbstractDofHandler) where {dim,T}
    grid = dh.grid

    _ndofs   = ndofs(part.element)
    _nnodes  = Ferrite.nnodes(getcelltype(part.element))
    nthreads = Threads.nthreads()
    nchunks = nthreads*10 #TODO: How many chunks?

    resize!(part.cache, nchunks)
    for i in 1:nchunks
        part.cache[i] = PartCache{dim,T}(_ndofs, _nnodes, part.element)
    end

    #TODO: BezierGrid does not work with create_incidence_matrix, so call them seperatly
    #threadsets = Ferrite.create_coloring(grid, part.cellset; ColoringAlgorithm.WorkStream)
    incidence_matrix = hot_fix_create_incidence_matrix(grid, part.cellset)
    threadsets = Ferrite.workstream_coloring(incidence_matrix, part.cellset)

    copy!(part.threadsets, threadsets)
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
    partstate::PartState,
    state::StateVariables{T},
    assemtype::ASSEMBLETYPE) where {dim,T,ET,MT}

    nchunks = length(part.cache)
    assemblers = [start_assemble(state.system_arrays.Kⁱ, state.system_arrays.fⁱ, fillzero=false) for _ in 1:nchunks]
   
    for color in part.threadsets
        Threads.@threads for (chunkrange, ichunk) in ChunkSplitters.chunks(color, nchunks) 
            for i in chunkrange
                cellid = color[i]
                lcellid = searchsortedfirst(part.cellset, cellid)
                
                partcache = part.cache[ichunk]
                assembler = assemblers[ichunk]

                _assemble_cell(part, partstate, partcache, cellid, lcellid,  assembler, state, dh, assemtype)
            end
        end
    end
            
    
end

function _assemble_cell(part::Part, partstate::PartState, cache::PartCache, cellid, lcellid, assembler, state, dh, assemtype)
    (; fe, ke, ue, due, Δue, coords, celldofs, element) = cache

    materialstate = partstate.materialstates[lcellid]
    cellstate     = partstate.elementstate[lcellid]
    stresses      = partstate.stresses[lcellid]
    strains       = partstate.strains[lcellid]

    fill!(fe, 0.0)
    (assemtype == STIFFMAT) && fill!(ke, 0.0)

    Ferrite.cellcoords!(coords, dh, cellid)
    Ferrite.celldofs!(celldofs, dh, cellid)

    Δue .= state.v[celldofs] 
    ue  .= state.d[celldofs]
    due .= state.v[celldofs]
    
    if assemtype == STIFFMAT
        integrate_forcevector_and_stiffnessmatrix!(element, cellstate, part.material, materialstate, stresses, strains, ke, fe, coords, Δue, ue, due, state.Δt)
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

function get_part_vtk_grid(part::Part)
    if part.geometry === nothing
        return nothing
    end
    return vtk_grid("mypart$(minimum(part.cellset))", part.geometry)
end

function eval_part_field_data(part::Part, dh, state, field_name::Symbol)
    fh = Five.getfieldhandler(dh, first(part.cellset))
    fieldidx = Ferrite.find_field(fh, field_name)
    if fieldidx === nothing #the field does not exist in this part
        return nothing 
    end

    field_offset = Ferrite.field_offset(fh, field_name) #TODO: change to fieldidx
    field_dim    = fh.fields[fieldidx].dim 
    space_dim = field_dim == 2 ? 3 : field_dim

    nnodes = getnnodes(dh.grid)
    data   = fill(Float64(NaN), space_dim, nnodes)
    Ferrite.reshape_field_data!(data, dh, state.d, field_offset, field_dim, part.cellset)
    return data[:, part.geometry.nodemapper]
end

function eval_part_node_data(part::Part, nodeoutput::VTKNodeOutput{<:StressOutput}, state, globaldata)
    #Extract stresses to interpolate
    qpdata = Vector{SymmetricTensor{2,3,Float64,6}}[]
    for (ic, cellid) in enumerate(part.cellset)
        stresses = state.partstates[cellid].stresses
        push!(qpdata, stresses)
    end
    data = zeros(SymmetricTensor{2,3,Float64,6}, getnnodes(globaldata.grid))
    _collect_nodedata!(data, part, qpdata, globaldata)
    return data[part.geometry.nodemapper]
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

function collect_nodedata!(data::Vector{FT}, part::Part{dim}, output::MaterialStateOutput{FT}, state::StateVariables{T}, globaldata) where {dim,FT,T}

    #Check if field exist in materialstate
    _cellid = first(part.cellset)
    first_state = first(state.partstates[_cellid].materialstates)
    if !hasproperty(first_state, output.field) 
        return
    end

    #Extract field to interpolate
    qpdata = Vector{FT}[] #TODO: allocate
    for (ic, cellid) in enumerate(part.cellset)
        matstates = state.partstates[cellid].materialstates
        field_states = getproperty.(matstates, output.field)
        push!(qpdata, field_states)
    end

    _collect_nodedata!(data, part, qpdata, globaldata)
end

function collect_nodedata!(data::Vector{FT}, part::Part{dim}, output::StressOutput, state::StateVariables{T}, globaldata) where {dim,FT,T}

    #Extract stresses to interpolate
    qpdata = Vector{SymmetricTensor{2,3,T,6}}[]
    for (ic, cellid) in enumerate(part.cellset)
        stresses = state.partstates[cellid].stresses
        push!(qpdata, stresses)
    end

    _collect_nodedata!(data, part, qpdata, globaldata)
end

function L2ProjectorByPassIGA(
    func_ip::Interpolation,
    grid::Ferrite.AbstractGrid;
    qr_lhs::QuadratureRule = Ferrite._mass_qr(func_ip),
    set = 1:getncells(grid),
    geom_ip::Interpolation = Ferrite.default_interpolation(typeof(grid.cells[first(set)])),
    qr_rhs::Union{QuadratureRule,Nothing}=nothing, # deprecated
)

    Ferrite._check_same_celltype(grid, collect(set)) # TODO this does the right thing, but gives the wrong error message if it fails

    fe_values_mass = CellScalarValues(qr_lhs, func_ip, geom_ip)

    # Create an internal scalar valued field. This is enough since the projection is done on a component basis, hence a scalar field.
    dh = MixedDofHandler(grid)
    field = Field(:_, func_ip, 1) # we need to create the field, but the interpolation is not used here
    fh = FieldHandler([field], Set(set))
    add!(dh, fh)
    _, vertex_dict, _, _ = Ferrite.__close!(dh)

    M = Ferrite._assemble_L2_matrix(fe_values_mass, set, dh)  # the "mass" matrix
    M_cholesky = cholesky(M)

    # For deprecated API
    fe_values = qr_rhs === nothing ? nothing :
                CellScalarValues(qr_rhs, func_ip, geom_ip)

    return L2Projector(func_ip, geom_ip, M_cholesky, dh, collect(set), vertex_dict[1], fe_values, qr_rhs)
end

function _collect_nodedata!(data::Vector{T}, part::Part{dim}, qpdata::Vector{Vector{FT}}, globaldata) where {dim,T,FT}
    
    (; grid, ) = globaldata

    #Set up quadrature rule
    celltype = getcelltype(part.element)
    geom_ip = Ferrite.default_interpolation(celltype)
    qr = getquadraturerule(part.element)

    projector = L2ProjectorByPassIGA(geom_ip, grid; set = part.cellset)
    projecteddata = project(projector, qpdata, qr; project_to_nodes=true); 

    #Reorder to the parts vtk
    for (ic, cellid) in enumerate(part.cellset)
        for nodeid in globaldata.grid.cells[cellid].nodes
            data[nodeid] = projecteddata[nodeid]
        end
    end
end

function collect_celldata!(data::Vector{FT}, part::Part{dim}, output::MaterialStateOutput, state::StateVariables{T}, globaldata) where {dim,FT,T}

    #Check if field exist in materialstate
    _cellid = first(part.cellset)
    first_state = first(state.partstates[_cellid].materialstates)
    if !hasproperty(first_state, output.field) 
        return
    end

    #Collect material state
    for (ic, cellid) in enumerate(part.cellset)
        materialstates = state.partstates[cellid].materialstates
        statevalues = getproperty.(materialstates, output.field)
        data[cellid] = output.func(statevalues)
    end
end

function collect_celldata!(data::Vector{FT}, part::Part{dim}, output::StressOutput, state::StateVariables{T}, globaldata) where {dim,FT, T}
    #Collect material state
    for (ic, cellid) in enumerate(part.cellset)
        stresses = state.partstates[cellid].stresses
        data[cellid] = output.func(stresses)
    end
end


function get_vtk_celldata(dh::Ferrite.AbstractDofHandler, part::Part, state::StateVariables) 
    asdf
    return nothing, nothing
end

function get_vtk_nodedata(dh::Ferrite.AbstractDofHandler, part::Part, state::StateVariables) 
    asdf
    return nothing, nothing
end
