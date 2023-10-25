export ExternalForceHandler, PointForce

abstract type AbstractExternalForce end

struct ExternalForceHandler{T,DH<:DofHandler}
    dh::DH
    external_forces::Vector{AbstractExternalForce}
end

function ExternalForceHandler{T}(dh::DofHandler) where T
    return ExternalForceHandler{T, typeof(dh)}(dh, AbstractExternalForce[])
end

function Ferrite.close!(ef::ExternalForceHandler)
    for (i, e) in enumerate(ef.external_forces)
        ForceType = typeof(e)
        ef.external_forces[i] = init_external_force!(e, ef.dh)
    end
end

function apply_external_forces!(dh, efh::ExternalForceHandler, state::StateVariables, globaldata)
    for force in efh.external_forces
        apply_external_force!(force, state, globaldata)
    end
end

"""
Point force
"""
struct PointForce <: AbstractExternalForce
    field::Symbol
    comps::Vector{Int}
    set::Vector{VertexIndex}
    func::Function
end

function PointForce(; field::Symbol, comps::AbstractVector{Int}, set, func::Function)
    return PointForce(field, collect(comps), collect(set), func)
end

function init_external_force!(force::PointForce, dh::DofHandler)
    Ferrite._check_same_celltype(dh.grid, cellid.(force.set))
    fh = getsubdofhandler(dh, cellid(first(force.set)))
    return PointForce(force.field, force.comps, force.set, force.func)
end

function apply_external_force!(force::PointForce, state::StateVariables, globaldata)
    for vertex in force.set
        sdh = getsubdofhandler(globaldata.dh, cellid(vertex))
        dofs = dofs_on_vertex(globaldata.dh, sdh, vertex, force.field, force.comps)
        X = Vec((0.0,0.0)) # TODO: position
        state.system_arrays.fᵉ[dofs] .= force.func(X, state.t)
    end
end


#
#
#
struct TractionForce{FV<:Ferrite.AbstractFaceValues} <: AbstractExternalForce
    set::Union{Vector{FaceIndex}, Vector{VertexIndex}, Vector{EdgeIndex}}
    traction::Function
    facevalues::FV
    udofs::UnitRange{Int}

    function TractionForce(set, traction::Function)
        return new{FaceValues}(set, traction)
    end

    function TractionForce(set, traction::Function, fv::FaceValues, udofs::UnitRange)
        return new{typeof(fv)}(set, traction, fv, udofs)
    end
end

function TractionForce(;
        set,
        traction::Function)

    return TractionForce(collect(set), traction)
end

function init_external_force!(force::TractionForce, dh::DofHandler)

    Ferrite._check_same_celltype(dh.grid, force.set)

    cellid, faceid = first(force.set)
    celltype = getcelltype(dh.grid, cellid)
    fh = getfieldhandler(dh, cellid)
    udofs = Ferrite.dof_range(fh, :u)
    i = Ferrite.find_field(fh, :u)
    i === nothing && error("The cellset containing cellid $(first(force.set)) does not have the field :u")
    
    field_ip = fh.fields[i].interpolation
    geo_ip = Ferrite.default_interpolation(celltype)
    face_ip = Ferrite.getlowerdim(field_ip)
    qr = Ferrite._mass_qr(face_ip)

    facevalues = FaceVectorValues(qr, field_ip, geo_ip)
    return TractionForce(force.set, force.traction, facevalues, udofs)
end

function apply_external_force!(ef::TractionForce{FV}, state::StateVariables{T}, globaldata) where {T,FV<:Ferrite.AbstractFaceValues}
    
    dh = globaldata.dh
    fv = ef.facevalues
    dim = Ferrite.getdim(dh)

    #Cache this?
    ndofs_field = getnbasefunctions(fv)
    _fcellid, fidx = first(ef.set)
    ndofs_cell = ndofs_per_cell(dh, _fcellid)
    ncoords = Ferrite.getngeobasefunctions(fv)

    fe = zeros(T, length(ef.udofs))
    coords = zeros(Vec{dim,T}, ncoords)
    celldofs = zeros(Int, ndofs_cell)
    total_area = 0.0
    for face in ef.set
        fill!(fe, 0.0)

        #extract cell id and face id 
        cellid, faceid = face

        #Get coords and dofs of cell
        Ferrite.cellcoords!(coords, dh, cellid)
        Ferrite.celldofs!(celldofs, dh, cellid)

        #Integrate
        area = _compute_external_traction_force!(fe, fv, coords, faceid, ef.traction, state.t)

        #Scatter
        total_area += area
        for (i,j) in enumerate(ef.udofs)
            state.system_arrays.fᵉ[celldofs[i]] += fe[j]
        end
    end
end


function _compute_external_traction_force!(fe::AbstractVector, fv::Ferrite.AbstractFaceValues, coords, faceid, traction::Function, time::T) where {T}

    reinit!(fv, coords, faceid)
    dA = 0
    for q_point in 1:getnquadpoints(fv)
        dΓ = getdetJdV(fv, q_point)

        X = spatial_coordinate(fv, q_point, coords)
        t = traction(X, time)

        dA += dΓ
        for i in 1:getnbasefunctions(fv)
            δui = shape_value(fv, q_point, i)
           
            fe[i] += (δui ⋅ t) * dΓ
        end
    end
   
    return dA
end
