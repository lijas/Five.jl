function IGALinearSolidElement(;
    thickness = 1.0, 
    qr_order::Int=2, 
    celltype::Type{<:IGA.BezierCell}, 
    dimstate::AbstractDim{dim} = MaterialModels.Dim{3}()
) where {dim}

    geom_ip = Ferrite.default_interpolation(celltype)
    ip = geom_ip
    
    qr = QuadratureRule{dim, RefCube}(qr_order)

    cv = IGA.BezierCellValues(CellVectorValues(qr, ip, geom_ip))
    return Five.LinearSolidElement{dim, (2,2,2), RefCube, Float64, typeof(cv), typeof(dimstate)}(thickness, celltype, cv, Field(:u, ip, dim), dimstate)
end

#Name it Part instead of Fe-part because it is the standard...
struct IGAPart{dim, T, E<:Five.AbstractElement, M<:MaterialModels.AbstractMaterial} <: Five.AbstractPart{dim}
    material::M
    cellset::Vector{Int}
    threadsets::Vector{Vector{Int}}
    element::E

    cache::Vector{Five.PartCache{dim,T,E}} #One for each thread
    geometry::Five.Optional{Any} # TODO: Where should the visualization geometry live?
end

function IGAPart(; 
    material::M,
    element::E,
    cellset,
    geometry=nothing
    ) where {E,M}

    dim = Ferrite.getdim(element)
    T = Float64
    _set = collect(cellset)
    sort!(_set) # YOLO
    return IGAPart{dim,T,E,M}(
        material, 
        _set,
        Vector{Int}[],
        element, 
        Five.PartCache{dim,T}[],
        geometry)
end

get_fields(part::IGAPart) = get_fields(part.element)
get_cellset(part::IGAPart) = part.cellset

function Five.construct_partstates(part::IGAPart{dim,T,ET,MT}) where {dim,T,ET,MT}

    ncells = length(part.cellset)
    nqp = Five.getnquadpoints(part.element)

    MaterialStateType = typeof( MaterialModels.initial_material_state(part.material) )
    ElementStateType = Five.elementstate_type(ET)

    states = Vector{Five.PartState{ElementStateType,MaterialStateType}}(undef, ncells)

    for i in 1:ncells
        #@show MaterialStateType,ElementStateType
        _materialstates = Vector{MaterialStateType}(undef, nqp)
        _elementstates = Vector{ElementStateType}(undef, nqp)
        for j in 1:nqp
            _materialstates[j] = MaterialModels.initial_material_state(part.material)
            _elementstates[j] = Five.initial_element_state(part.element)
        end

        states[i] = Five.PartState(_elementstates, _materialstates, zeros(SymmetricTensor{2,3,Float64,6}, nqp), zeros(SymmetricTensor{2,3,Float64,6}, nqp))
    end
    return states
end

function hot_fix_create_incidence_matrix(g::Ferrite.AbstractGrid, cellset)
    cell_containing_node = Dict{Int, Set{Int}}()
    for cellid in cellset
        cell = Ferrite.getcells(g, cellid)
        for v in cell.nodes
            _set = get!(Set{Int}, cell_containing_node, v)
            push!(_set, cellid)
        end
    end

    I, J, V = Int[], Int[], Bool[]
    for (_, cells) in cell_containing_node
        for cell1 in cells # All these cells have a neighboring node
            for cell2 in cells
                # if true # cell1 != cell2
                if cell1 != cell2
                    push!(I, cell1)
                    push!(J, cell2)
                    push!(V, true)
                end
            end
        end
    end

    incidence_matrix = sparse(I, J, V, getncells(g), getncells(g))
    return incidence_matrix
end

function Five.init_part!(part::IGAPart{dim, T}, dh::Ferrite.AbstractDofHandler) where {dim,T}
    grid = dh.grid

    _ndofs = ndofs(part.element)
    _nnodes = Ferrite.nnodes(Ferrite.getcelltype(part.element))

    nthreads = Threads.nthreads()
    
    resize!(part.cache, nthreads)
    for i in 1:nthreads
        part.cache[i] = Five.PartCache{dim,T}(_ndofs, _nnodes, part.element)
    end

    #Hot fix cells for parts with only one cell (bar_example.jl)
    #TODO: this is fixed in latest ferrite version
    local threadsets
    if length(part.cellset) == 1
        threadsets = Vector{Int}[[first(part.cellset)]]
    else
        #TODO: BezierGrid does not work with create_incidence_matrix, so call them seperatly
        #threadsets = Ferrite.create_coloring(grid, part.cellset; ColoringAlgorithm.WorkStream)
        incidence_matrix = Five.hot_fix_create_incidence_matrix(grid, part.cellset)
        threadsets = Ferrite.workstream_coloring(incidence_matrix, part.cellset)
    end

    copy!(part.threadsets, threadsets)
end

function Five.assemble_stiffnessmatrix_and_forcevector!(dh::Ferrite.AbstractDofHandler, 
    part::IGAPart,
    state::Five.StateVariables)

    Five._assemble_part!(dh, part,state, Five.STIFFMAT)

end

function Five.assemble_forcevector!(dh::Ferrite.AbstractDofHandler, 
    part::Part,
    state::StateVariables)

    Five. _assemble_part!(dh, part,state, Five.FORCEVEC)

end

function Five.assemble_fstar!(dh::Ferrite.AbstractDofHandler, 
    part::Part,
    state::StateVariables)

    Five._assemble_part!(dh, part,state, Five.FSTAR)

end

function Five.assemble_dissipation!(
    dh    ::Ferrite.AbstractDofHandler, 
    part  ::Part,
    state ::StateVariables)

    if !(Five.is_dissipative(part.material) || Five.is_dissipative(part.element))
        return 
    end

    Five._assemble_part!(dh, part, state, Five.DISSI)

end

function Five._assemble_part!(dh::Ferrite.AbstractDofHandler, 
    part::IGAPart{dim,T,ET,MT},
    state::StateVariables{T},
    assemtype::Five.ASSEMBLETYPE) where {dim,T,ET,MT}

    assemblers = [start_assemble(state.system_arrays.Kⁱ, state.system_arrays.fⁱ, fillzero=false) for _ in 1:Threads.nthreads()]
   
    ElementState = Five.elementstate_type(ET)
    MaterialState = Five.materialstate_type(MT)

    Δt = state.Δt

    coordscahce = [IGA.getcoordinates(dh.grid, first(first(part.threadsets))) for _ in 1:Threads.nthreads()]
    for tset in part.threadsets
        Threads.@threads :static for cellid in tset

            cache = part.cache[Threads.threadid()]
            assembler = assemblers[Threads.threadid()]
            coords = coordscahce[Threads.threadid()]
            
            (; fe, ke, ue, due, Δue, celldofs, element) = cache
        
            partstate::PartState{ElementState, MaterialState} = state.partstates[cellid]

            materialstate = partstate.materialstates
            cellstate     = partstate.elementstate
            stresses      = partstate.stresses
            strains       = partstate.strains

            fill!(fe, 0.0)
            (assemtype == Five.STIFFMAT) && fill!(ke, 0.0)

            Ferrite.getcoordinates!(coords, dh.grid, cellid)

            Ferrite.celldofs!(celldofs, dh, cellid)

            Δue .= state.v[celldofs] #Dont need both Δue and (v and Δt)
            ue .= state.d[celldofs]
            due .= state.v[celldofs]
            
            if assemtype == Five.STIFFMAT
                Five.integrate_forcevector_and_stiffnessmatrix!(element, cellstate, part.material, materialstate, stresses, strains, ke, fe, coords, Δue, ue, due, Δt)
                assemble!(assembler, celldofs, fe, ke)
            elseif assemtype == Five.FORCEVEC
                Five.integrate_forcevector!(element, cellstate, part.material, materialstate, fe, coords, Δue, ue, due, Δt)
                state.system_arrays.fⁱ[celldofs] += fe
            elseif assemtype == Five.FSTAR
                error("Broken code, fix")
                prev_partstate::Five.get_partstate_type(part) = state.prev_partstates[cellid]
                prev_materialstate = prev_partstate.materialstates

                Five.integrate_fstar!(element, cellstate, part.material, prev_materialstate, fe, coords, Δue, ue, due, Δt)
                state.system_arrays.fᴬ[celldofs] += fe
            elseif assemtype == Five.DISSI
                ge = Base.RefValue(zero(T))
                Five.integrate_dissipation!(element, cellstate, part.material, materialstate, fe, ge, coords, Δue, ue, due, Δt)
                state.system_arrays.fᴬ[celldofs] += fe
                state.system_arrays.G[] += ge[]
            end

        end
    end
    
end



function Five.get_part_vtk_grid(part::IGAPart)
    if part.geometry === nothing
        return nothing
    end
    return vtk_grid("mypart$(minimum(part.cellset))", part.geometry)
end

function Five.eval_part_field_data(part::IGAPart, dh, state, field_name::Symbol)
    data = _evaluate_at_geometry_nodes!(part.geometry, dh, state.d, field_name, part.cellset)
    return data
end

Ferrite._mass_qr(::IGA.BernsteinBasis{3, (2, 2, 2)}) = QuadratureRule{3, RefCube}(3)
function Five.eval_part_node_data(part::IGAPart, nodeoutput::Five.VTKNodeOutput{<:Five.StressOutput}, state, globaldata)
    #Extract stresses to interpolate
    qpdata = Vector{SymmetricTensor{2,3,Float64,6}}[]
    for (ic, cellid) in enumerate(part.cellset)
        stresses = state.partstates[cellid].stresses
        for i in 1:length(stresses)
            stresses[i] = zero(SymmetricTensor{2,3,Float64,6})
        end
        push!(qpdata, stresses)
    end
    
    #
    #
    #Set up quadrature rule
    celltype = getcelltype(part.element)
    geom_ip = Ferrite.default_interpolation(celltype)
    qr = getquadraturerule(part.element)
    vtktype = Ferrite.cell_to_vtkcell(celltype)
    reorder = IGA.Ferrite_to_vtk_order(celltype)

    projector = L2Projector(geom_ip, globaldata.grid; set = part.cellset)
    projecteddata = IGA.igaproject(projector, qpdata, qr; project_to_nodes=true); 
    #
    #


    #
    #
    local_node_coords = Ferrite.reference_coordinates(geom_ip)
    n_eval_points = length(local_node_coords)

    offset = 0
    data = zeros(6, size(part.geometry.coords,2))
    for cellid in part.geometry.cellset

        cellnodes = globaldata.grid.cells[cellid].nodes
        nodevalues = projecteddata[collect(cellnodes)]
        IGA._distribute_vtk_point_data!(globaldata.grid.beo[cellid], data, nodevalues, offset)
        #vtk_cellnodes = collect((1:n_eval_points) .+ offset)
        #for (nodeid, vtknodeid) in zip(cellnodes, vtk_cellnodes)
        #    σ = projecteddata[nodeid]
        #    data[1:6, vtknodeid] .= tovoigt(σ)
        #end

        offset += n_eval_points
    end

    #
    #
    #

    return data
end

function eval_part_node_data(part::IGAPart, nodeoutput::VTKNodeOutput{MaterialStateOutput{MaterialState_t}}, state, globaldata) where MaterialState_t
    
    #Extract stresses to interpolate
    _cellid = first(part.cellset)
    first_state = first(state.partstates[_cellid].materialstates)
    if !hasproperty(first_state, nodeoutput.type.field) 
        return
    end
    
    qpdata = Vector{MaterialState_t}[]
    for (ic, cellid) in enumerate(part.cellset)
        matstates = state.partstates[cellid].materialstates
        field_states = getproperty.(matstates, nodeoutput.type.field)
        push!(qpdata, field_states)
    end
    
    #
    #
    #Set up quadrature rule
    celltype = getcelltype(part.element)
    geom_ip = Ferrite.default_interpolation(celltype)
    qr = getquadraturerule(part.element)
    vtktype = Ferrite.cell_to_vtkcell(celltype)
    reorder = IGA.Ferrite_to_vtk_order(celltype)

    projector = L2Projector(geom_ip, globaldata.grid; set = part.cellset)
    projecteddata = IGA.igaproject(projector, qpdata, qr; project_to_nodes=true); 
    #
    #


    #
    #
    local_node_coords = Ferrite.reference_coordinates(geom_ip)
    n_eval_points = length(local_node_coords)

    offset = 0
    ncomp = length(first(projecteddata))
    data = zeros(ncomp, size(part.geometry.coords,2))
    for cellid in part.geometry.cellset

        cellnodes = globaldata.grid.cells[cellid].nodes
        nodevalues = projecteddata[collect(cellnodes)]
        IGA._distribute_vtk_point_data!(globaldata.grid.beo[cellid], data, nodevalues, offset)
        #vtk_cellnodes = collect((1:n_eval_points) .+ offset)
        #for (nodeid, vtknodeid) in zip(cellnodes, vtk_cellnodes)
        #    σ = projecteddata[nodeid]
        #    data[1:6, vtknodeid] .= tovoigt(σ)
        #end

        offset += n_eval_points
    end

    #
    #
    #

    return data
end

function Five.post_part!(dh, part::IGAPart, states::StateVariables)
    
end

function Five.commit_part!(dh::Ferrite.AbstractDofHandler, part::IGAPart, state::StateVariables)
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

function Five.get_vtk_field(dh::Ferrite.AbstractDofHandler, part::IGAPart{dim,T}, state::StateVariables, field_name::Symbol) where {dim,T}
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

function Five._get_vtk_field!(data::Matrix, dh::Ferrite.AbstractDofHandler, part::IGAPart{dim,T}, state::StateVariables, offset::Int, nvars::Int) where {dim,T}

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

function Five.collect_nodedata!(data::Vector{FT}, part::IGAPart{dim}, output::MaterialStateOutput{FT}, state::StateVariables{T}, globaldata) where {dim,FT,T}

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

function Five.collect_nodedata!(data::Vector{FT}, part::IGAPart{dim}, output::Five.StressOutput, state::StateVariables{T}, globaldata) where {dim,FT,T}

    #Extract stresses to interpolate
    qpdata = Vector{SymmetricTensor{2,3,T,6}}[]
    for (ic, cellid) in enumerate(part.cellset)
        stresses = state.partstates[cellid].stresses
        push!(qpdata, stresses)
    end

    _collect_nodedata!(data, part, qpdata, globaldata)
end

function _collect_nodedata!(data::Vector{T}, part::IGAPart{dim}, qpdata::Vector{Vector{FT}}, globaldata) where {dim,T,FT}
    
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

function Five.collect_celldata!(data::Vector{FT}, part::IGAPart{dim}, output::MaterialStateOutput, state::StateVariables{T}, globaldata) where {dim,FT,T}

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

function Five.collect_celldata!(data::Vector{FT}, part::IGAPart{dim}, output::Five.StressOutput, state::StateVariables{T}, globaldata) where {dim,FT, T}
    #Collect material state
    for (ic, cellid) in enumerate(part.cellset)
        stresses = state.partstates[cellid].stresses
        data[cellid] = output.func(stresses)
    end
end


function Five.get_vtk_celldata(dh::Ferrite.AbstractDofHandler, part::IGAPart, state::StateVariables) 
    asdf
    return nothing, nothing
end

function Five.get_vtk_nodedata(dh::Ferrite.AbstractDofHandler, part::IGAPart, state::StateVariables) 
    asdf
    return nothing, nothing
end




#
struct IGASubGridGeometry
    grid::Ferrite.AbstractGrid
    cellset::Vector{Int}
    vtkcells::Vector{WriteVTK.MeshCell}
    coords::Matrix{Float64}
end 

function IGASubGridGeometry(grid::Ferrite.AbstractGrid{dim}, cellset) where dim

    #
    vtkcells = MeshCell[]
	coords_vec = Vec{dim,Float64}[]
    weights = Float64[]
	cellorders = Int[]
	offset = 0
    celltype = typeof(grid.cells[first(cellset)])
    vtktype = Ferrite.cell_to_vtkcell(celltype)
    reorder = IGA.Ferrite_to_vtk_order(celltype)

    n = getnbasefunctions(Ferrite.default_interpolation(celltype))
    xb = zeros(Vec{dim}, n)
    wb = zeros(Float64, n)
    x  = zeros(Vec{dim}, n)
    w  = zeros(Float64, n)

    cellset = collect(cellset) #Collect cellset becuase order matters

    for cellid in cellset
		cell = grid.cells[cellid]

		IGA.get_bezier_coordinates!(xb, wb, x, w, grid, cellid)

	    append!(coords_vec, xb)
		append!(weights, wb)

		cellnodes = (1:length(cell.nodes)) .+ offset
		offset += length(cell.nodes)

        push!(vtkcells, MeshCell(vtktype, collect(cellnodes[reorder])))
    end
    
    coords_mat = reshape(reinterpret(Float64, coords_vec), (dim, length(coords_vec)))

    return IGASubGridGeometry(grid, cellset, vtkcells, coords_mat)
end

function WriteVTK.vtk_grid(filename, geometry::IGASubGridGeometry)
    return vtk_grid(filename, geometry.coords, geometry.vtkcells)
end


function _evaluate_at_geometry_nodes!(
    geometry  ::IGASubGridGeometry, 
    dh        ::Ferrite.AbstractDofHandler, 
    a         ::Vector, 
    fieldname ::Symbol, 
    compcells ::Union{Set{Int},Vector{Int}})

    fieldname ∈ Ferrite.getfieldnames(dh) || error("Field $fieldname not found in the dofhandler.")
    field_dim = Ferrite.getfielddim(dh, fieldname)
    space_dim = field_dim == 2 ? 3 : field_dim
    
    nviznodes = size(geometry.coords, 2)
    data = fill(Float64(NaN), space_dim, nviznodes)  # set default value

    # Assume that only one of the field handlers can "project" to the geometry
    local fh
    for _fh in dh.fieldhandlers
        if first(compcells) in _fh.cellset
            fh = _fh
            @assert Set{Int}(compcells) == _fh.cellset
            break
        end
    end

    # Check if this fh contains this field, otherwise continue to the next
    field_idx = Ferrite.find_field(fh, fieldname)
    field_idx === nothing && error("The field handler does not contain the field $fieldname")
    ip     = Ferrite.getfieldinterpolation(fh, field_idx)
    drange = Ferrite.dof_range(fh, fieldname)

    _evaluate_at_geometry_nodes!(data, dh, a, ip, drange, field_dim, geometry)
        
    return data
end


function _evaluate_at_geometry_nodes!(data, dh, a, ip, drange, field_dim, geometry::IGASubGridGeometry)

    local_node_coords = Ferrite.reference_coordinates(ip)
    n_eval_points = length(local_node_coords)
    qr = QuadratureRule{field_dim,RefCube}(zeros(length(local_node_coords)), local_node_coords)
    cv = IGA.BezierCellValues( CellVectorValues(qr, ip) )

    offset = 0
    for cellid in geometry.cellset

        coords = getcoordinates(dh.grid, cellid)
        reinit!(cv, coords)

        dofs = zeros(Int, ndofs_per_cell(dh, cellid))
        celldofs!(dofs, dh, cellid)
        ae = a[dofs[drange]]

        cellnodes = (1:n_eval_points) .+ offset

        for iqp in 1:n_eval_points
            u = function_value(cv, iqp, ae)
            data[1:field_dim, cellnodes[iqp]] .= u
        end

        offset += n_eval_points
    end

    return data
end

