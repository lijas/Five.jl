export ExternalForceHandler, PointForce

abstract type AbstractExternalForce end

struct ExternalForceHandler{T}
    external_forces::Vector{AbstractExternalForce}
end

function ExternalForceHandler{T}() where T
    return ExternalForceHandler{T}(AbstractExternalForce[])
end

function JuAFEM.close!(ef::ExternalForceHandler, dh::MixedDofHandler)
    for (i, e) in enumerate(ef.external_forces)
        ForceType = typeof(e)
        ef.external_forces[i] = init_external_force!(e, dh)
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
    JuAFEM._check_same_celltype(dh.grid, cellid.(force.set))
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
struct TractionForce{FV<:JuAFEM.Values} <: AbstractExternalForce
    field::Symbol
    comps::Vector{Int}
    set::Union{Vector{FaceIndex}, Vector{VertexIndex}, Vector{EdgeIndex}}
    traction::Function

    facevalues::FV
end

function TractionForce{dim}(;
    set, 
    traction::Function) where {dim}

    return TractionForce{dim,T,FV}(collect(faces), traction, facevalues, beo)
end

function init_external_force!(force::TractionForce, dh::MixedDofHandler)
    
    id = cellid(first(force.set))

    return TractionForce(force.field, force.comps, force.set, force.func, fh)
end

function apply_external_force!(dh::JuAFEM.AbstractDofHandler, ef::TractionForce{FV}, state::StateVariables, globaldata) where {FV<:JuAFEM.Values}
    
    fv = ef.facevalues
    ndofs = getnbasefunctions(fv)
    fe = zeros(T, ndofs)
    ncoords = JuAFEM.getngeobasefunctions(fv)
    coords = zeros(Vec{dim,T}, ncoords)
    celldofs = zeros(Int, ndofs)

    total_area = 0
    for (i, face) in enumerate(ef.faces)
        fill!(fe, 0.0)

        #extract cell id and face id 
        cellid,faceid = (face[1], face[2])
        #Get coords and dofs of cell
        JuAFEM.cellcoords!(coords, dh, cellid)
        JuAFEM.celldofs!(celldofs, dh, cellid)
        
        #Special case isogeomtric parts
        if FV <: IGA.BezierValues
            Ce = ef.bezier_operators[i]
            IGA.set_bezier_operator!(fv, Ce)
            coords .= IGA.compute_bezier_points(Ce, coords)
        end
        
        area = _compute_external_traction_force!(fv, coords, faceid, ef.traction, fe, state.t)
        total_area +=area
        system_arrays.fᵉ[celldofs] += fe
    end
    @show total_area
end


function _compute_external_traction_force!(fv::JuAFEM.Values{dim,T}, cellcoords, faceid, traction::Function, fe::AbstractVector, t::T) where {dim,T}

    reinit!(fv, cellcoords, faceid)
    dA = 0
    for q_point in 1:getnquadpoints(fv)
        dΓ = getdetJdV(fv, q_point)

        X = spatial_coordinate(fv, q_point, cellcoords)
        t =traction(X,time)

        dA += dΓ
        for i in 1:getnbasefunctions(fv)
            δui = shape_value(fv, q_point, i)
           
            fe[i] += (δui ⋅ t) * dΓ
        end
    end
   
    return dA
end
