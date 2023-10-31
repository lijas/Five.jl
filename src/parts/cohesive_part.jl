
const CZPart{dim,T} = Part{dim,T,<:CohesiveElement}

#=
function eval_part_field_data(geometry::SubGridGeometry, part::CZPart, dh, state, fieldname::Symbol)
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
    @show ip
    qr = QuadratureRule{shape}(zeros(length(local_node_coords)), local_node_coords)
    cv = SurfaceVectorValues(qr, ip.ip)

    drange = dof_range(sdh, field_idx)

    # Function barrier
    ##
    ##
    ue = zeros(Float64, length(drange))
    if RT <: Vec && cv isa CellValues{<:ScalarInterpolation}
        uer = reinterpret(RT, ue)
    else
        uer = ue
    end
    for cell in CellIterator(sdh)
        # Note: We are only using the shape functions: no reinit!(cv, cell) necessary
        @assert getnquadpoints(cv) == length(cell.nodes)
        for (i, I) in pairs(drange)
            ue[i] = state.dofs[cell.dofs[I]]
        end
        for (qp, nodeid) in pairs(cell.nodes)
            val = mid_surf_value(cv, qp, uer)
            if data isa Matrix # VTK
                data[1:length(val), nodeid] .= val
                data[(length(val)+1):end, nodeid] .= 0 # purge the NaN
            else
                data[nodeid] = val
            end
        end
    end

    return data[:, geometry.nodemapper]
end=#


function collect_nodedata!(data::Vector{FT}, part::Part{dim,T,<:CohesiveElement}, output::StressOutput, state::StateVariables{T}, globaldata) where {dim,FT,T} 
end

function default_geometry(part::CZPart, grid)
    return nothing#SubGridGeometry(grid, part.cellset)
end
