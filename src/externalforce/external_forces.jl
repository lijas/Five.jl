export ExternalForceHandler, PointForce

abstract type AbstractExternalForce end

struct ExternalForceHandler{T,DH<:Ferrite.MixedDofHandler}
    dh::DH
    external_forces::Vector{AbstractExternalForce}
end

function ExternalForceHandler{T}(dh::MixedDofHandler) where T
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
    
    fieldhandler::FieldHandler
end

function PointForce(; field::Symbol, comps::AbstractVector{Int}, set, func::Function)
    return PointForce(field, collect(comps), collect(set), func, FieldHandler())
end

function init_external_force!(force::PointForce, dh::MixedDofHandler)
    Ferrite._check_same_celltype(dh.grid, cellid.(force.set))
    fh = getfieldhandler(dh, cellid(first(force.set)))
    return PointForce(force.field, force.comps, force.set, force.func, fh)
end

function apply_external_force!(force::PointForce, state::StateVariables, globaldata)
    for vertex in force.set
        dofs = dofs_on_vertex(globaldata.dh, force.fieldhandler, vertex, force.field, force.comps)
        X = Vec((0.0,0.0)) # TODO: position
        state.system_arrays.fᵉ[dofs] .= force.func(X, state.t)
    end
end


#
#
#
struct TractionForce{FV<:Ferrite.Values} <: AbstractExternalForce
    field::Symbol
    set::Union{Vector{FaceIndex}, Vector{VertexIndex}, Vector{EdgeIndex}}
    traction::Function

    facevalues::FV
end

function TractionForce(;
        field::Symbol,
        set,
        traction::Function,
        celltype::Type{<:Ferrite.AbstractCell{dim}}) where dim

    ip = Ferrite.default_interpolation(celltype)
    fip = Ferrite.getlowerdim(ip)
    qr = Ferrite._mass_qr(fip)
    fv = FaceVectorValues(qr, ip)

    val = traction(zero(Vec{dim}), 0.0)
    @assert length(val) == 3

    return TractionForce(field, collect(set), traction, fv)
end

function init_external_force!(force::TractionForce, dh::MixedDofHandler)
    return force
end

function apply_external_force!(ef::TractionForce{FV}, state::StateVariables{T}, globaldata) where {T,FV<:Ferrite.Values}
    
    dh = globaldata.dh
    fv = ef.facevalues
    dim = Ferrite.getdim(dh)

    #Cache this?
    ndofs = getnbasefunctions(fv)
    fe = zeros(T, ndofs)
    ncoords = Ferrite.getngeobasefunctions(fv)
    coords = zeros(Vec{dim,T}, ncoords)
    celldofs = zeros(Int, ndofs)

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
        state.system_arrays.fᵉ[celldofs] += fe
    end
end


function _compute_external_traction_force!(fe::AbstractVector, fv::Ferrite.Values{dim,T}, coords, faceid, traction::Function, t::T) where {dim,T}

    reinit!(fv, coords, faceid)
    dA = 0
    for q_point in 1:getnquadpoints(fv)
        dΓ = getdetJdV(fv, q_point)

        X = spatial_coordinate(fv, q_point, coords)
        t = traction(X, t)

        dA += dΓ
        for i in 1:getnbasefunctions(fv)
            δui = shape_value(fv, q_point, i)
           
            fe[i] += (δui ⋅ t) * dΓ
        end
    end
   
    return dA
end
