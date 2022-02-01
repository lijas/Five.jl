export SolidElement
"""
SolidElement{dim,order,shape,T,CV<:Ferrite.Values}

"""

struct SolidElement{dim,order,shape,T,CV<:Ferrite.Values} <: AbstractElement
    thickness::T #used in 2d

    celltype::Type{<:Cell}
    cv::CV
    field::Field

    dimension::AbstractDim
end

Ferrite.getnquadpoints(e::SolidElement) = getnquadpoints(e.cv)
Ferrite.ndofs(e::SolidElement) = getnbasefunctions(e.cv)
Ferrite.getcelltype(e::SolidElement) = e.celltype
getncoords(s::SolidElement) = Ferrite.getngeobasefunctions(s.cv)

has_constant_massmatrix(::SolidElement) = true

get_fields(e::SolidElement{dim,order,shape,T}) where {dim,order,shape,T} = return [e.field]

function SolidElement{dim,order,refshape,T}(;thickness = 1.0, qr_order::Int=2, celltype::Type{<:Cell}, dimension = ThreeD()) where {dim, order, refshape, T}
    
    ip = Lagrange{dim, refshape, order}()
    geom_ip = Ferrite.default_interpolation(celltype)

    qr = QuadratureRule{dim, refshape}(qr_order)

    cv = CellVectorValues(qr, ip, geom_ip)
    return SolidElement{dim,order,refshape,T,typeof(cv)}(thickness, celltype, cv, Field(:u, ip, dim), dimension)
end

function calculate_minimum_timestep(element::SolidElement{2,1,RefCube,T,M}, material::AbstractMaterial, cell::CellIterator, ue::Vector, due::Vector) where {T,M}
    
    L1 = norm(cell.coords[1] - cell.coords[2])
    L2 = norm(cell.coords[2] - cell.coords[3])
    L3 = norm(cell.coords[3] - cell.coords[4])
    L4 = norm(cell.coords[4] - cell.coords[1])

    Area = 0.5*(cell.coords[1][1]*cell.coords[2][2] + 
                cell.coords[2][1]*cell.coords[3][2] + 
                cell.coords[3][1]*cell.coords[4][2] + 
                cell.coords[4][1]*cell.coords[1][2] - 
                cell.coords[2][1]*cell.coords[1][2] - 
                cell.coords[3][1]*cell.coords[2][2] -
                cell.coords[4][1]*cell.coords[3][2] - 
                cell.coords[1][1]*cell.coords[4][2])

    characteristic_length = Area = max(L1,L2,L3,L4)

    density = density(material)
    stiffness = linear_stiffness(material.E)

    wave_speed = sqrt(stiffness/density)
    min_timestep = characteristic_length/wave_speed
    return min_timestep 
end

function integrate_fstar!(element::SolidElement{dim, order, shape, T}, 
    elementstate::AbstractElementState, 
    material::AbstractMaterial, 
    materialstate::AbstractVector{<:AbstractMaterialState},
    fe::Vector, 
    cell, 
    Δue::Vector,
    ue::Vector,
    due::Vector,
    Δt::T) where {dim, order, shape, T}
    
    integrate_forcevector_and_stiffnessmatrix!(element, elementstate,
     material, materialstate, zeros(T, length(fe), length(fe)), fe, cell, Δue, ue, due, Δt)
end

function integrate_forcevector!(element::SolidElement{dim, order, shape, T}, 
        elementstate::AbstractElementState, 
        material::AbstractMaterial, 
        materialstate::AbstractVector{<:AbstractMaterialState},
        fe::Vector, 
        cell, 
        Δue::Vector,
        ue::Vector,
        due::Vector,
        Δt::T) where {dim, order, shape, T}


        cv = element.cv
        reinit!(cv, cell)
        ndofs = Ferrite.ndofs(element)
    
        δE = zeros(SymmetricTensor{2, dim, eltype(ue)}, ndofs)
    
        for qp in 1:getnquadpoints(cv)
            ∇u = function_gradient(cv, qp, ue)
            dΩ = getdetJdV(cv, qp) * element.thickness
    
            # strain and stress + tangent
            F = one(∇u) + ∇u
            E = symmetric(1/2 * (F' ⋅ F - one(F)))
    
            S, _, new_matstate = material_response(material.dimension, material, F, materialstate[qp])
            materialstate[qp] = new_matstate

            for i in 1:ndofs
                δFi = shape_gradient(cv, qp, i)
                δE[i] = symmetric(1/2*(δFi'⋅F + F'⋅δFi))
            end
    
            for i in 1:ndofs
                fe[i] += (δE[i] ⊡ S) * dΩ
            end
        end
end

function integrate_massmatrix!(element::SolidElement{dim, order, shape, T, M}, elstate::AbstractElementState, material::AbstractMaterial, cell, me::Matrix, ue::AbstractVector, due::AbstractVector) where {dim, order, shape, T, M}

    cv = element.cv
    reinit!(cv, cell)
    ndofs = Ferrite.ndofs(element)

    for qp in 1:getnquadpoints(cv)
        dV = getdetJdV(cv, qp) * element.thickness
        for i in 1:ndofs
            Ni = shape_value(cv, qp, i)
            for j in 1:ndofs
                Nj = shape_value(cv, qp, j)
                me[i, j] += density(material)*Ni⋅Nj*dV;
            end
        end
    end
end

function integrate_forcevector_and_stiffnessmatrix!(element::SolidElement{dim, order, shape, T}, 
    elementstate::AbstractElementState, 
    material::AbstractMaterial, 
    materialstate::AbstractVector{<:AbstractMaterialState},
    ke::AbstractMatrix, 
    fe::Vector, 
    cell, 
    Δue::Vector,
    ue::Vector,
    due::Vector,
    Δt::T) where {dim, order, shape, T}

    cv = element.cv
    reinit!(cv, cell)
    ndofs = Ferrite.ndofs(element)

    δE = zeros(SymmetricTensor{2, dim, eltype(ue)}, ndofs)

    for qp in 1:getnquadpoints(cv)
        ∇u = function_gradient(cv, qp, ue)
        dΩ = getdetJdV(cv, qp) * element.thickness

        # strain and stress + tangent
        F = one(∇u) + ∇u
        E = symmetric(1/2 * (F' ⋅ F - one(F)))

        #assume output is S and ∂S∂E for now
        C = tdot(F)
        S, dSdE, new_matstate = material_response(element.dimension, ∂S∂E(), material, F, materialstate[qp])
        materialstate[qp] = new_matstate
        #U = sqrt(C)
        #R = F⋅inv(U)

        #∂S∂E = otimesu(R,R) ⊡ _∂S∂E ⊡ otimesu(R',R')
        #S = R ⋅ _S ⋅R'

        # Hoist computations of δE
        for i in 1:ndofs
            δFi = shape_gradient(cv, qp, i)
            δE[i] = symmetric(1/2*(δFi'⋅F + F'⋅δFi))
        end

        for i in 1:ndofs
            δFi = shape_gradient(cv, qp, i)
            δu = shape_value(cv, qp, i)
            fe[i] += (δE[i] ⊡ S) * dΩ
            δE∂S∂E = δE[i] ⊡ dSdE
            S∇δu = S ⋅ δFi'
            for j in 1:ndofs
                δ∇uj = shape_gradient(cv, qp, j)
                ke[i, j] += (δE∂S∂E ⊡ δE[j] + S∇δu ⊡ δ∇uj' ) * dΩ
            end
        end
    end
end


function integrate_dissipation!(
    element       :: SolidElement{dim_p,dim_s,CV},
    elementstate  :: AbstractElementState,
    material      :: AbstractMaterial,
    materialstate :: AbstractArray{<:AbstractMaterialState},
    fe            :: Vector{T},
    ge            :: Base.RefValue{T},
    coords, 
    Δue           :: Vector,
    ue            :: Vector,
    due           :: Vector,
    Δt            :: T
    ) where {dim_p,dim_s,CV,T}
    

    #if !is_dissipative(material)
    #    return
    #end

    cv = element.cv
    reinit!(cv, coords)
    ndofs = Ferrite.ndofs(element)

    for qp in 1:getnquadpoints(cv)

        dΩ = getdetJdV(cv, qp) * element.thickness

        ∇u = function_gradient(cv, qp, ue)
        F = one(∇u) + ∇u
        C = tdot(F)
        g, ∂g∂C = constitutive_driver_dissipation(material, C, materialstate[qp])

        ge[] += g * dΩ

        for i in 1:ndofs
            δFi = shape_gradient(cv, qp, i)
            δCi = δFi'⋅F + F'⋅δFi
            fe[i] += (∂g∂C ⊡ δCi) * dΩ
        end
    end
end
