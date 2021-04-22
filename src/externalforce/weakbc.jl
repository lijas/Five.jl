#
#
#
struct WeakBoundaryCondition{dim,T,FV<:JuAFEM.Values} <: AbstractExternalForce
    faces::IndexVectors
    prescribed_displacement::Function
    local_dofs::Vector{Int}
    facevalues::FV
end

function WeakBoundaryCondition{dim,T}(dh::JuAFEM.AbstractDofHandler, faces, prescribed_displacement::Function, components, facevalues::FV) where {dim,T,FV<:JuAFEM.Values{dim,T}}
   # @assert(length(prescribed_displacement(0.0,0.0))==length(components))

    celltype = typeof(dh.grid.cells[first(faces)[1]])

    for cell in dh.grid.cells[2:end]
        if typeof(cell) != celltype
            error("Cells must be of same type")
        end 
    end

    local_dofs = Int[]
    for i in 1:getnbasefunctions(facevalues)÷dim
        for d in 1:dim
            if d in components
                push!(local_dofs, (i-1)*dim+d)
            end
        end
    end

    return WeakBoundaryCondition{dim,T,FV}(faces, prescribed_displacement, local_dofs, facevalues)
end

function apply_external_force!(dh::JuAFEM.AbstractDofHandler, ef::WeakBoundaryCondition{dim,T,FV}, state::StateVariables, prev_state::StateVariables, system_arrays::SystemArrays, globaldata) where {dim,T,FV<:JuAFEM.Values}
    
    fv = ef.facevalues
    local_dofs = ef.local_dofs

    ndofs = getnbasefunctions(fv)
    fe = zeros(T, ndofs)
    ke = zeros(T, ndofs, ndofs)
    ncoords = getnbasefunctions(fv)÷dim
    coords = zeros(Vec{dim,T}, ncoords)
    celldofs = zeros(Int, ndofs)

    _g = ef.prescribed_displacement(zero(Vec{dim,T}), state.t)
    g = Vec((_g...,))
    
    A = 0.0
    for (i, face) in enumerate(ef.faces)
        fill!(fe, 0.0)
        fill!(ke, 0.0)

        #extract cell id and face id 
        cellid,faceid = (face[1], face[2])
        
        #Get coords and dofs of cell
        cellcoords!(coords, dh, cellid)
        JuAFEM.celldofs!(celldofs, dh, cellid)
        lcd = celldofs[local_dofs]
        
        A += _compute_weak_boundary_condition!(fv, coords, faceid, g, ke, fe, state.d[celldofs])

        system_arrays.fᵉ[lcd] += fe[local_dofs]
        system_arrays.Kᵉ[lcd,lcd] += ke[local_dofs,local_dofs]
    end
    @show A

end


function _compute_weak_boundary_condition!(fv::JuAFEM.Values{dim,T}, coords::AbstractVector{Vec{dim,T}}, faceid, prescr_disp::Function, ke, fe::AbstractVector, ue::AbstractVector) where {dim,T}
    reinit!(fv, cellcoords, faceid)
    stiffness = 1e7
    dA = 0.0
    for q_point in 1:getnquadpoints(fv)
        dΓ = getdetJdV(fv, q_point)
        u = function_value(fv, q_point, ue)
        
        X = spatial_coordinate(fv, q_point, coords)
        g = prescr_disp(X,time)

        dA += dΓ
        for i in 1:getnbasefunctions(fv)
            δui = shape_value(fv, q_point, i)  
            
            fe[i] += stiffness*(u-g)⋅δui * dΓ
            for j in 1:getnbasefunctions(fv)
                δuj = shape_value(fv, q_point, j)
                ke[i,j] += stiffness*δui⋅δuj * dΓ
            end
        end
    end
    return dA

end