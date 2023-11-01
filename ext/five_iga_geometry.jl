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
    reorder = IGA._iga_to_vtkorder(celltype)

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

function Five.eval_part_field_data(geometry::IGASubGridGeometry, part::Part, dh, state, fieldname::Symbol)
    data = _evaluate_at_geometry_nodes!(geometry, dh, state.d, fieldname, part.cellset)
    return data
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
    for _fh in dh.subdofhandlers
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

    refshape = Ferrite.getrefshape(ip)
    @show typeof(ip)
    local_node_coords = Ferrite.reference_coordinates(ip)
    n_eval_points = length(local_node_coords)
    qr = QuadratureRule{refshape}(zeros(length(local_node_coords)), local_node_coords)
    cv = CellValues( qr, ip, ip.ip)

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
            data[(field_dim+1):end, cellnodes[iqp]] .= 0.0 # purge the NaN
        end

        offset += n_eval_points
    end

    return data
end


#Ferrite._mass_qr(::IGA.BernsteinBasis{3, (2, 2, 2)}) = QuadratureRule{3, RefCube}(3)
#Ferrite._mass_qr(::IGA.BernsteinBasis{3, (2, 2, 2)}) = QuadratureRule{3, RefCube}(3)
function Five.eval_part_node_data(geometry::IGASubGridGeometry, part::IGAPart, partstate::Five.PartState, nodeoutput::Five.VTKNodeOutput{<:Five.StressOutput}, state, globaldata)
    #Extract stresses to interpolate
    qpdata = Vector{SymmetricTensor{2,3,Float64,6}}[]
    for (ic, cellid) in enumerate(part.cellset)
        stresses = partstate.stresses[cellid]
        #for i in 1:length(stresses)
        #    stresses[i] = zero(SymmetricTensor{2,3,Float64,6})
        #end
        push!(qpdata, stresses)
    end
    #
    #
    #Set up quadrature rule
    celltype = getcelltype(part.element)
    geom_ip = Ferrite.default_interpolation(celltype)
    qr = Five.getquadraturerule(part.element)
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

function Five.eval_part_node_data(part::IGAPart, nodeoutput::VTKNodeOutput{MaterialStateOutput{MaterialState_t}}, state, globaldata) where MaterialState_t
    
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
