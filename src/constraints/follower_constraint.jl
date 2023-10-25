

struct FollowerConstraint
    faces::Set{FaceIndex}
    mastervertex::VertexIndex
    field_name::Symbol
    component::Int
end

#Keyword constructor
function FollowerConstraint(;
    faces::Set{FaceIndex},
    mastervertex::VertexIndex,
    field::Symbol,
    component::Int,)

    return FollowerConstraint(faces, mastervertex, field, component)
end


function Ferrite.add!(ch::ConstraintHandler, c::FollowerConstraint)
    cellid, faceid = first(c.faces) #assume all cells in the cellset are the same
    fh            = getsubdofhandler(ch.dh, cellid)
    field_idx     = Ferrite.find_field(fh, c.field_name)
    interpolation = Ferrite.getfieldinterpolation(fh, field_idx)
    field_dim     = Ferrite.getfielddim(fh, field_idx)
    offset        = Ferrite.field_offset(fh, c.field_name)

    _add!(ch, c, interpolation, field_dim, offset)
    return ch
end

function _add!(ch::ConstraintHandler, c::FollowerConstraint, interpolation::Interpolation, field_dim::Int, offset::Int)

    dh = ch.dh
    component = c.component

    #
    local_face_dofs  , local_face_dofs_offset   = Ferrite._local_face_dofs_for_bc(interpolation, field_dim, [component], offset)
    local_vertex_dofs, local_vertex_dofs_offset = Ferrite._local_face_dofs_for_bc(interpolation, field_dim, [component], offset, Ferrite.vertices)
    
    #
    ncelldofs = ndofs_per_cell(dh, Ferrite.cellid(first(c.faces)))
    cdofs = zeros(Int, ncelldofs)
    dofs_on_face = Int[]
    for (cellid, faceidx) in c.faces
        celldofs!(cdofs, dh, cellid)
        r = local_face_dofs_offset[faceidx]:(local_face_dofs_offset[faceidx+1]-1)
        append!(dofs_on_face, cdofs[local_face_dofs[r]]) 
    end
    
    unique!(dofs_on_face)

    
    #Get master dof
    cellid, vertexidx = c.mastervertex 
    celldofs!(cdofs, dh, cellid)
    r = local_vertex_dofs_offset[vertexidx]#:(local_vertex_dofs_offset[vertexidx+1]-1)
    @assert (local_vertex_dofs[r] |> length) == 1
    masterdof = cdofs[local_vertex_dofs[r]][1]

    filter!(i->i!==masterdof, dofs_on_face)
    
    for dof in dofs_on_face
        ac = Ferrite.AffineConstraint(dof, [masterdof => 1.0], 0.0)
        add!(ch, ac)
    end

end