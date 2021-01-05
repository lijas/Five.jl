

function IGAPart{dim,T}(material::AbstractMaterial, partset::Vector{Int}, element::AbstractElement, Ce::AbstractVector) where {dim,T}

    plot_orders = ntuple(_->1, dim)
    
    points = JuAFEM.reference_coordinates(IGA.BernsteinBasis{dim,plot_orders}())

    weights = zeros(T, length(points))
    qr = QuadratureRule{dim,RefCube,T}(weights,points)
    
    ip = element.fields[1].interpolation #Assume
    cv_plot = CellScalarValues(qr,ip)
    
    return IGAPart{dim,T,typeof(element),typeof(material),}(material, partset, element, Ce, cv_plot)
end

get_bezier_operator(part::IGAPart, ie::Int) = return part.Cb[ie]

function init_part!(part::IGAPart, dh::JuAFEM.AbstractDofHandler) 

end

#Could be combined with "FE"-part assemble-loop because they are very similar
#=function assemble_stiffnessmatrix_and_forcevector!(dh::JuAFEM.AbstractDofHandler, 
                                                        part::IGAPart{dim,T}, 
                                                        state::StateVariables,
                                                        prev_state::StateVariables,
                                                        system_arrays::SystemArrays) where {dim,T}
    assembler = start_assemble(system_arrays.Kⁱ, system_arrays.fⁱ, fillzero=false)
    element = part.element

    fe = zeros(T, ndofs(element))
    ke = zeros(T, ndofs(element), ndofs(element))
    coords = zeros(Vec{dim,T}, getncoords(element))
    bezier_coords = zeros(Vec{dim,T}, getncoords(element))
    bezier_ue = zeros(Vec{dim,T}, getncoords(element))
    celldofs = zeros(Int, ndofs(element))

    for (localid, cellid) in enumerate(part.cellset)

        cellcoords!(coords, dh, cellid)
        JuAFEM.celldofs!(celldofs, dh, cellid)

        #Update bezier shit
        Ce = get_bezier_operator(part, localid)
        IGA.set_bezier_operator!(part.element.cv, Ce)
        bezier_coords .= IGA.compute_bezier_points(Ce, coords)

        new_materialstates =      state.materialstates[cellid]
        materialstate      = prev_state.materialstates[cellid]
        cellstate = prev_state.cellstates[cellid]
        
        fill!(fe, 0.0)
        fill!(ke, 0.0)
        
        Δue = state.d[celldofs] - prev_state.d[celldofs]
        ue = state.d[celldofs]
        due = state.v[celldofs]
        
        dV   = integrate_forcevector_and_stiffnessmatrix!(element, cellstate, part.elementinfo, part.material, materialstate, new_materialstates, ke, fe, bezier_coords, Δue, ue, due)
        
        assemble!(assembler, celldofs, fe, ke)
    end

end=#



function get_vtk_grid(dh::JuAFEM.AbstractDofHandler, part::IGAPart{dim,T}) where {dim,T}
    celltype = getcelltype(part.element)

    n_plot_points_dims = ntuple(_->2, dim)

    cv = part.cv_plot #Use the element CellValues. In future might store its own in IGApart instead

    cls = MeshCell[]
    node_coords = Vec{dim,T}[]
    node_offset = 0
    
    for (local_id,celldata) in enumerate(CellIterator2(dh, part.element, part.cellset))
        cellid = celldata.current_cellid[]

        Ce = part.Cb[local_id]
        bezier_coords = IGA.compute_bezier_points(Ce, celldata.coords)

        for point_id in 1:getnquadpoints(cv)
            X = spatial_coordinate(cv, point_id, bezier_coords)
            push!(node_coords, X)
        end

        n_plot_points = prod(n_plot_points_dims)::Int
        nodeind2nodeid = reshape(1:n_plot_points, n_plot_points_dims)
        indeces = CartesianIndices(n_plot_points_dims .- 1 )[:]
        addons  = CartesianIndices(Tuple(fill(2,dim)))[:]

        for index in indeces

            nodeids = Int[]
            for a in 1:2^dim
                nodendex = Tuple(index).+ Tuple(addons[a]) .-1
                _nodeid = nodeind2nodeid[nodendex...]
                push!(nodeids, _nodeid)
            end
            
            VTK_CELL = (dim==2) ? JuAFEM.VTKCellTypes.VTK_QUAD : JuAFEM.VTKCellTypes.VTK_HEXAHEDRON
            nodeids = (dim==2) ? nodeids[[1,2,4,3]] : nodeids[[1,2,4,3,5,6,8,7]]
            
            push!(cls, MeshCell(VTK_CELL, nodeids .+ node_offset))

        end 
        node_offset = length(node_coords)
    end
   
    return cls, node_coords
end

function get_vtk_displacements(dh::JuAFEM.AbstractDofHandler, part::IGAPart{dim,T}, state::StateVariables) where {dim,T}
    celltype = getcelltype(part.element)

    #Only accept elements with displacement field (:u)
    @assert(length(part.element.fields) == 1 && part.element.fields[1].name == :u)

    node_coords = Vec{dim,T}[]
    cv = part.cv_plot
    bcv = IGA.BezierValues(cv)

    for (local_id,celldata) in enumerate(CellIterator2(dh, part.element, part.cellset))
        ue = state.d[celldata.celldofs]
        ue_vec = reinterpret(Vec{dim,T}, ue)

        cellid = celldata.current_cellid[]

        Ce = get_bezier_operator(part, local_id)        
        ue_bezier = IGA.compute_bezier_points(Ce, ue_vec)
        
        IGA.set_bezier_operator!(bcv, Ce)

        bezier_coords = IGA.compute_bezier_points(Ce, celldata.coords)
        reinit!(bcv, bezier_coords)

        for point_id in 1:getnquadpoints(cv)
            u2 = function_value(cv, point_id, ue_bezier)
            u = function_value(bcv, point_id, ue_vec)
            @assert(u≈u2)
            
            push!(node_coords, u)
        end

    end
    
    return node_coords
end