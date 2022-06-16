"""
SolidElement{dim,order,shape,T,CV<:Ferrite.Values}
"""

struct SolidElement{
        dim,
        order,
        shape,
        T          <: AbstractFloat,
        CV         <: Ferrite.Values,
        DIM       <: MaterialModels.AbstractDim
    } <: AbstractElement

    thickness::T #used in 2d

    celltype::Type{<:Cell}
    cv::CV
    field::Field
    dimstate::DIM
end

Ferrite.getnquadpoints(e::SolidElement) = getnquadpoints(e.cv)
Ferrite.ndofs(e::SolidElement) = getnbasefunctions(e.cv)
Ferrite.getcelltype(e::SolidElement) = e.celltype
getncoords(s::SolidElement) = Ferrite.getngeobasefunctions(s.cv)
has_constant_massmatrix(::SolidElement) = true
get_fields(e::SolidElement)= return [e.field]

function SolidElement{dim,order,refshape,T}(;
        thickness::Float64 = 1.0, 
        qr_order::Int=2, 
        celltype::Type{<:Cell}, #Todo: Should be obtained from the grid instead...
        dimstate::AbstractDim{dim} = MaterialModels.Dim{3}()) where {dim, order, refshape, T}
    
    ip = Lagrange{dim, refshape, order}()
    geom_ip = Ferrite.default_interpolation(celltype)

    qr = QuadratureRule{dim, refshape}(qr_order)

    cv = CellVectorValues(qr, ip, geom_ip)
    return SolidElement{dim, order, refshape, T, typeof(cv), typeof(dimstate)}(thickness, celltype, cv, Field(:u, ip, dim), dimstate)
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

        S, ∂S∂E, new_matstate = material_response(element.dimstate, material, E, materialstate[qp])
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

function integrate_massmatrix!(
    element::SolidElement{dim, order, shape, T, M}, 
    elstate::AbstractElementState, 
    material::AbstractMaterial, 
    coords, 
    me::Matrix, 
    ue::AbstractVector, 
    due::AbstractVector) where {dim, order, shape, T, M}

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
        
        #@assert( strainmeasure(material) == MaterialModels.GreenLagrange() )
        S, ∂S∂E, new_matstate = material_response(element.dimstate, material, E, materialstate[qp])
        materialstate[qp] = new_matstate

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