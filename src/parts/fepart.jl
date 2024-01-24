
struct PartCache{dim,T,E,coords_t}
    ue::Vector{T}
    due::Vector{T}
    Δue::Vector{T}
    fe::Vector{T}
    ke::Matrix{T}
    celldofs::Vector{Int}
    
    #Make the coords a type parameter for IGA
    #For standard fe-problems, coords::Vector{Vec{dim,T}}
    coords::coords_t #Vec{dim,T}
    element::E
end

function create_part_cache(element::E) where E<:AbstractElement
    dim = Ferrite.getdim(element)
    T = Float64
    _ndofs   = ndofs(element)
    _nnodes  = Ferrite.nnodes(getcelltype(element))

    coords_t = Vector{Vec{dim,T}}
    PartCache{dim,T,E,coords_t}(
        zeros(T,_ndofs), 
        zeros(T,_ndofs), 
        zeros(T,_ndofs), 
        zeros(T,_ndofs), 
        zeros(T,_ndofs,_ndofs), 
        zeros(Int,_ndofs), 
        zeros(Vec{dim,T},_nnodes),
        deepcopy(element)
    )
end

#Name it Part instead of Fe-part because it is the standard...
struct Part{dim, T, E<:AbstractElement, M<:AbstractMaterial, cache_t<:PartCache{dim,T,E}} <: AbstractPart{dim}
    material::M
    cellset::Vector{Int}
    threadsets::Vector{Vector{Int}}
    element::E

    cache::Vector{cache_t} #One for each thread
end

function Part(; 
    material::M,
    element::E,
    cellset,
    ) where {E,M}

    dim = Ferrite.getdim(element)
    T = Float64

    _set = collect(cellset)
    sort!(_set) # YOLO

    coords_t = Vector{Vec{dim,T}}
    cache_t = PartCache{dim,T,E,coords_t}

    return Part{dim,T,E,M,cache_t}(
        material, 
        _set,
        Vector{Int}[],
        element, 
        cache_t[])
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
    Ferrite._check_same_celltype(grid, part.cellset)

    nthreads = Threads.nthreads()
    nchunks = nthreads*10 #TODO: How many chunks?

    resize!(part.cache, nchunks)
    for i in 1:nchunks
        part.cache[i] = create_part_cache(part.element)
    end

    #incidence_matrix = hot_fix_create_incidence_matrix(grid, part.cellset)
    #Ferrite.workstream_coloring(incidence_matrix, part.cellset)
    threadsets = Ferrite.create_coloring(grid, part.cellset; alg=ColoringAlgorithm.WorkStream) 
    copy!(part.threadsets, threadsets)
end

@enum ASSEMBLETYPE FORCEVEC STIFFMAT FSTAR DISSI

function assemble_stiffnessmatrix_and_forcevector!(dh::Ferrite.AbstractDofHandler, 
    part::Part,
    partstate::PartState,
    state::StateVariables)

    _assemble_part!(dh, part, partstate, state, STIFFMAT)

end

function assemble_forcevector!(dh::Ferrite.AbstractDofHandler, 
    part::Part,
    partstate::PartState,
    state::StateVariables)

    _assemble_part!(dh, part, partstate, state, FORCEVEC)

end

function assemble_fstar!(dh::Ferrite.AbstractDofHandler, 
    part::Part,
    partstate::PartState,
    state::StateVariables)

    _assemble_part!(dh, part, partstate, state, FSTAR)

end

function assemble_dissipation!(
    dh    ::Ferrite.AbstractDofHandler, 
    part::Part,
    partstate::PartState,
    state ::StateVariables)

    if !(is_dissipative(part.material) || is_dissipative(part.element))
        return 
    end

    _assemble_part!(dh, part, partstate, state, DISSI)

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

function _assemble_cell(part::Part, partstate::PartState, cache::PartCache, cellid, lcellid, assembler, state::StateVariables{T}, dh, assemtype) where T
    (; fe, ke, ue, due, Δue, coords, celldofs, element) = cache

    materialstate = partstate.materialstates[lcellid]
    cellstate     = partstate.elementstate[lcellid]
    stresses      = partstate.stresses[lcellid]
    strains       = partstate.strains[lcellid]

    fill!(fe, 0.0)
    (assemtype == STIFFMAT) && fill!(ke, 0.0)

    Ferrite.getcoordinates!(coords, dh.grid, cellid)
    Ferrite.celldofs!(celldofs, dh, cellid)

    Δue .= state.v[celldofs] 
    ue  .= state.d[celldofs]
    due .= state.v[celldofs]
    
    if assemtype == STIFFMAT
        integrate_forcevector_and_stiffnessmatrix!(element, cellstate, part.material, materialstate, stresses, strains, ke, fe, coords, Δue, ue, due, state.Δt)
        assemble!(assembler, celldofs, fe, ke)
    elseif assemtype == FORCEVEC
        integrate_forcevector!(element, cellstate, part.material, materialstate, fe, coords, Δue, ue, due, state.Δt)
        state.system_arrays.fⁱ[celldofs] += fe
    elseif assemtype == FSTAR
        error("Broken code, fix")
        #integrate_fstar!(element, cellstate, part.material, prev_materialstate, fe, coords, Δue, ue, due, state.Δt)
        #state.system_arrays.fᴬ[celldofs] += fe
    elseif assemtype == DISSI
        ge = Base.RefValue(zero(T))
        integrate_dissipation!(element, cellstate, part.material, materialstate, fe, ge, coords, Δue, ue, due, state.Δt)
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

function default_geometry(part::Part, grid)
    return SubGridGeometry(grid, part.cellset)
end

function eval_part_field_data(geometry::SubGridGeometry, part::Part, dh, state, fieldname::Symbol)
    sdh = Five.getsubdofhandler(dh, first(part.cellset))

    #
    field_idx = Ferrite.find_field(dh, fieldname)
    ip = Ferrite.getfieldinterpolation(dh, field_idx)
    RT = ip isa ScalarInterpolation ? Float64 : Vec{Ferrite.n_components(ip),Float64}

    # VTK output of solution field (or L2 projected scalar data)
    n_c = Ferrite.n_components(ip)
    vtk_dim = n_c == 2 ? 3 : n_c # VTK wants vectors padded to 3D
    data = fill(NaN * zero(Float64), vtk_dim, getnnodes( Ferrite.get_grid(dh)))

    # Check if this sdh contains this field, otherwise continue to the next
    field_idx = Ferrite._find_field(sdh, fieldname)
    @assert field_idx !== nothing

    # Set up CellValues with the local node coords as quadrature points
    CT = Ferrite.getcelltype(sdh)
    ip = Ferrite.getfieldinterpolation(sdh, field_idx)
    ip_geo = Ferrite.default_interpolation(CT)
    local_node_coords = Ferrite.reference_coordinates(ip_geo)
    shape = Ferrite.getrefshape(ip)

    qr = QuadratureRule{shape}(zeros(length(local_node_coords)), local_node_coords)
    if ip isa VectorizedInterpolation
        # TODO: Remove this hack when embedding works...
        cv = CellValues(qr, ip.ip, ip_geo)
    else
        cv = CellValues(qr, ip, ip_geo)
    end
    drange = dof_range(sdh, field_idx)

    # Function barrier
    Ferrite._evaluate_at_grid_nodes!(data, sdh, state.d, cv, drange, RT)

    return data[:, geometry.nodemapper]
end

function eval_part_node_data(geometry::SubGridGeometry, part::Part, partstate::PartState, nodeoutput::VTKNodeOutput{<:StressOutput}, state, globaldata)
    #Extract stresses to interpolate
    qpdata = Vector{SymmetricTensor{2,3,Float64,6}}[]
    for lcellid in 1:length(part.cellset)
        stresses = partstate.stresses[lcellid]
        push!(qpdata, stresses)
    end
    data = zeros(SymmetricTensor{2,3,Float64,6}, getnnodes(globaldata.grid))
    _collect_nodedata!(data, part, qpdata, globaldata)
    return data[geometry.nodemapper]
end


function eval_part_cell_data(geometry::SubGridGeometry, part::Part, partstate::PartState, nodeoutput::VTKCellOutput{<:StressOutput}, state, globaldata)
    #Extract stresses to interpolate
    celldata = SymmetricTensor{2,3,Float64,6}[]
    for lcellid in 1:length(part.cellset)
        stresses = partstate.stresses[lcellid]
        meanstess = mean(stresses)
        push!(celldata, meanstess)
    end
    return celldata
end

function eval_part_cell_data(geometry::SubGridGeometry, part::Part, partstate::PartState, output::VTKCellOutput{MaterialStateOutput{DT}}, state, globaldata) where DT
    
    #Check if field exist in materialstate
    first_state = first(first(partstate.materialstates))
    if !hasproperty(first_state, output.type.field) 
        return
    end

    #Extract stresses to interpolate
    celldata = DT[]
    for (ic, cellid) in enumerate(part.cellset)
        matstates = partstate.materialstates[ic]
        field_states = getproperty.(matstates, output.type.field)
        mean_field = mean(field_states)
        push!(celldata, mean_field)
    end

    return celldata
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


function collect_nodedata!(data::Vector{FT}, part::Part{dim}, partstate::PartState, output::MaterialStateOutput{FT}, state::StateVariables{T}, globaldata) where {dim,FT,T}

    #Check if field exist in materialstate
    first_state = first(first(partstate.materialstates))
    if !hasproperty(first_state, output.field) 
        return
    end

    #Extract field to interpolate
    qpdata = Vector{FT}[] #TODO: allocate
    for (ic, cellid) in enumerate(part.cellset)
        matstates = partstate.materialstates[ic]
        field_states = getproperty.(matstates, output.field)
        push!(qpdata, field_states)
    end

    _collect_nodedata!(data, part, qpdata, globaldata)
end

function collect_nodedata!(data::Vector{FT}, part::Part{dim}, partstate::PartState, output::StressOutput, state::StateVariables{T}, globaldata) where {dim,FT,T}

    #Extract stresses to interpolate
    qpdata = Vector{SymmetricTensor{2,3,T,6}}[]
    for (ic, cellid) in enumerate(part.cellset)
        stresses = partstate.stresses[ic]
        push!(qpdata, stresses)
    end

    _collect_nodedata!(data, part, qpdata, globaldata)
end


function _collect_nodedata!(data::Vector{T}, part::Part{dim}, qpdata::Vector{Vector{FT}}, globaldata) where {dim,T,FT}
    
    (; grid, ) = globaldata

    #Set up quadrature rule
    celltype = getcelltype(part.element)
    geom_ip = Ferrite.default_interpolation(celltype)
    qr = getquadraturerule(part.element)

    projector = L2Projector(geom_ip, grid; set = part.cellset)
    projecteddata = project(projector, qpdata, qr); 

    projection_at_nodes = evaluate_at_grid_nodes(projector, projecteddata)
    data .= projection_at_nodes 
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
