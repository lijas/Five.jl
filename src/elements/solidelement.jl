"""
SolidElement
"""

struct SolidElement{sdim,CV<:Ferrite.AbstractCellValues, DIMSTATE<:MaterialModels.AbstractDim} <: AbstractElement{sdim}
    celltype::Type{<:Ferrite.AbstractCell}
    cv::CV
    dimstate::DIMSTATE
    thickness::Float64
end

initial_element_state(::SolidElement)  = EmptyElementState()

getquadraturerule(e::SolidElement) = _getquadraturerule(e.cv)
Ferrite.getnquadpoints(e::SolidElement) = getnquadpoints(e.cv)
Ferrite.ndofs(e::SolidElement) = getnbasefunctions(e.cv)
Ferrite.getcelltype(e::SolidElement) = e.celltype
has_constant_massmatrix(::SolidElement) = true
get_fields(e::SolidElement) = [(:u, e.cv.ip)]

function SolidElement(;
        celltype::Type{<:Ferrite.AbstractCell{refshape}}, 
        thickness                   = 1.0, 
        qr_order::Int               = 2, 
        dimstate::AbstractDim{sdim} = MaterialModels.Dim{3}(),
        ip::Interpolation           = Ferrite.default_interpolation(celltype),
    ) where {refshape, sdim}

    @assert Ferrite.getrefshape(ip) == refshape "refshape of element does not match that of the given interpolation"

    geo_ip = Ferrite.default_interpolation(celltype)
    qr = QuadratureRule{refshape}(qr_order)
    cv = CellValues(qr, ip^sdim, geo_ip^sdim)

    return SolidElement{sdim}(celltype, cv, dimstate, thickness)
end

function integrate_fstar!(element::SolidElement, 
    elementstate::Vector{<:AbstractElementState}, 
    material::AbstractMaterial, 
    materialstate::AbstractVector{<:AbstractMaterialState},
    fe::Vector, 
    cell, 
    Δue::Vector,
    ue::Vector,
    due::Vector,
    Δt::AbstractFloat) 
    
    integrate_forcevector_and_stiffnessmatrix!(element, elementstate,
     material, materialstate, zeros(T, length(fe), length(fe)), fe, cell, Δue, ue, due, Δt)
end

function integrate_forcevector!(element::SolidElement{dim}, 
        elementstate::Vector{<:AbstractElementState}, 
        material::AbstractMaterial, 
        materialstate::AbstractVector{<:AbstractMaterialState},
        stresses::Vector{<:SymmetricTensor{2,3,T}},
        strains::Vector{<:SymmetricTensor{2,3,T}},
        fe::Vector, 
        cell, 
        Δue::Vector,
        ue::Vector,
        due::Vector,
        Δt::T) where {dim, T}


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

        S, ∂S∂E, new_matstate = material_response(element.dimstate, material, E, materialstate[qp])
        materialstate[qp] = new_matstate

        #Store stress and strain (always as 3d)
        stresses[qp] = dim==3 ? S : MaterialModels.increase_dim(S)
        strains[qp] = dim==3 ? S : MaterialModels.increase_dim(E)

        for i in 1:ndofs
            δFi = shape_gradient(cv, qp, i)
            δE[i] = symmetric(1/2*(δFi'⋅F + F'⋅δFi))
        end

        for i in 1:ndofs
            fe[i] += (δE[i] ⊡ S) * dΩ
        end
    end
end

function integrate_massmatrix!(
    element::SolidElement{dim}, 
    elementstate::Vector{<:AbstractElementState}, 
    material::AbstractMaterial, 
    coords, 
    me::Matrix, 
    ue::AbstractVector, 
    due::AbstractVector) where {dim}

    cv = element.cv
    reinit!(cv, coords)
    ndofs = Ferrite.ndofs(element)

    for qp in 1:getnquadpoints(cv)
        dV = getdetJdV(cv, qp) * element.thickness
        for i in 1:ndofs
            Ni = shape_value(cv, qp, i)
            for j in 1:ndofs
                Nj = shape_value(cv, qp, j)
                me[i, j] += 1.0*Ni⋅Nj*dV;
            end
        end
    end
end

function integrate_forcevector_and_stiffnessmatrix!(element::SolidElement{dim}, 
    elementstate::Vector{<:AbstractElementState}, 
    material::AbstractMaterial, 
    materialstate::AbstractVector{<:AbstractMaterialState},
    stresses::Vector{<:SymmetricTensor{2,3,T}},
    strains::Vector{<:SymmetricTensor{2,3,T}},
    ke::AbstractMatrix, 
    fe::Vector, 
    cell, 
    Δue::Vector,
    ue::Vector,
    due::Vector,
    Δt::T) where {dim, T}

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
        
        #@assert( strainmeasure(material) == MaterialModels.GreenLagrange() )
        S, ∂S∂E, new_matstate = material_response(element.dimstate, material, E, materialstate[qp])
        materialstate[qp] = new_matstate

        #Store stress and strain (always as 3d)
        stresses[qp] = dim==3 ? S : MaterialModels.increase_dim(S)
        strains[qp] = dim==3 ? S : MaterialModels.increase_dim(E)

        # Hoist computations of δE
        for i in 1:ndofs
            δFi = shape_gradient(cv, qp, i)
            δE[i] = symmetric(1/2*(δFi'⋅F + F'⋅δFi))
        end

        for i in 1:ndofs
            δFi = shape_gradient(cv, qp, i)
            δu = shape_value(cv, qp, i)
            fe[i] += (δE[i] ⊡ S) * dΩ
            δE∂S∂E = δE[i] ⊡ ∂S∂E
            S∇δu = S ⋅ δFi'
            for j in 1:ndofs
                δ∇uj = shape_gradient(cv, qp, j)
                ke[i, j] += (δE∂S∂E ⊡ δE[j] + S∇δu ⊡ δ∇uj' ) * dΩ
            end
        end
    end
end


function integrate_dissipation!(
    element       :: SolidElement,
    elementstate  :: Vector{<:AbstractElementState}, 
    material      :: AbstractMaterial,
    materialstate :: AbstractArray{<:AbstractMaterialState},
    fe            :: Vector{T},
    ge            :: Base.RefValue{T},
    coords, 
    Δue           :: Vector,
    ue            :: Vector,
    due           :: Vector,
    Δt            :: AbstractFloat
    ) where T

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