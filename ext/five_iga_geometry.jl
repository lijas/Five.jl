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

function eval_part_field_data(geometry::IGASubGridGeometry, part::Part, dh, state, fieldname::Symbol)
    data = _evaluate_at_geometry_nodes!(part.geometry, dh, state.d, fieldname, part.cellset)
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
            data[(field_dim+1):end, cellnodes[iqp]] .= 0.0 # purge the NaN
        end

        offset += n_eval_points
    end

    return data
end

