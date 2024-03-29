export CohesiveElement

"""
CohesiveElement

Any order, any shape, any dim, 

"""

struct CohesiveElement{dim_s,CV} <: AbstractElement{dim_s}
    thickness2d::Float64

    field::Field
    celltype::Type{<:CohesiveCell}
    cv::CV
end

initial_element_state(::CohesiveElement) = EmptyElementState()

getquadraturerule(e::CohesiveElement) = e.cv.qr
Ferrite.getnquadpoints(e::CohesiveElement) = getnquadpoints(e.cv)
Ferrite.ndofs(e::CohesiveElement) = getnbasefunctions(e.cv)
Ferrite.getcelltype(e::CohesiveElement) = e.celltype
has_constant_massmatrix(::CohesiveElement) = true
get_fields(e::CohesiveElement) = return [e.field]


function CohesiveElement(;
        thickness::Float64 = 1.0, 
        order::Int, 
        nqp::Int = order+1,
        celltype::Type{<:CohesiveCell{dim_s}}
        ) where {dim_s}

    ip = CohesiveZoneInterpolation(Lagrange{dim_s-1,RefCube,order}())
    
    geom_ip = Ferrite.default_interpolation(celltype)
    mid_qr = QuadratureRule{dim_s-1,RefCube}(nqp)

    cv = SurfaceVectorValues(mid_qr, ip, geom_ip)

    return CohesiveElement{dim_s,typeof(cv)}(thickness, Field(:u, ip, dim_s), celltype, cv)
end

function integrate_forcevector_and_stiffnessmatrix!(
        element::CohesiveElement{dim_s,CV}, 
        elementstate::Vector{<:AbstractElementState}, 
        material::AbstractMaterial, 
        materialstate::AbstractArray{<:AbstractMaterialState}, 
        stresses::Vector{<:SymmetricTensor{2,3,T}},
        strains::Vector{<:SymmetricTensor{2,3,T}},
        ke::AbstractMatrix, 
        fe::Vector{T}, 
        cell, 
        Δue::Vector,
        ue::Vector,
        due::Vector,
        Δt::T) where {dim_s,CV,T}
    
    cv = element.cv
    
    ndofs = Ferrite.ndofs(element)

    xe = cell + reinterpret(Vec{dim_s,T}, ue)#[1:length(cell)*dim_s])

    reinit!(cv, xe)

    A = 0
    for qp in 1:getnquadpoints(cv)
        
        #Rotation matrix
        R = getR(cv,qp)
        @assert(det(R) ≈ 1.0)
        
        dΓ = getdetJdA(cv, qp) * element.thickness2d
        A += dΓ
        
        J = function_value(cv, qp, ue)
        Ĵ = R'⋅J
        
        #constitutive_driver
        t̂, ∂t∂Ĵ, new_matstate = material_response(material, Ĵ, materialstate[qp])
        materialstate[qp] = new_matstate
        #stresses[q_point] = t
        #strains[q_point] = J

        t = R ⋅ t̂
        ∂t∂J = R ⋅ ∂t∂Ĵ ⋅ R'

        for i in 1:ndofs
            δui = shape_value(cv, qp, i)

            fe[i] += (t ⋅ δui) * dΓ
            for j in 1:ndofs
                δuj = shape_value(cv, qp, j)
                ke[i,j] += δui⋅∂t∂J⋅δuj * dΓ
            end
        end
    end

end

function integrate_fstar!(element::CohesiveElement{dim_s,CV}, 
        elementstate::AbstractElementState, 
        material::AbstractMaterial, 
        materialstate::AbstractArray{<:AbstractMaterialState}, 
        fe::Vector{T}, 
        cell, 
        Δue::Vector,
        ue::Vector,
        due::Vector,
        Δt::T) where {dim_s,CV,T}
    
    error("Needs fixing")
    cv = element.cv
    
    ndofs = Ferrite.ndofs(element)

    xe = cell + reinterpret(Vec{dim_s,T}, ue)

    reinit!(cv, xe)

    A = 0
    for qp in 1:getnquadpoints(cv)
        
        #Rotation matrix
        R = getR(cv,qp)
        #@assert(det(R)==1.0)

        dΓ = getdetJdA(cv, qp) * element.thickness2d
        A += dΓ
        
        J = function_value(cv, qp, ue)
        Ĵ = R'⋅J
        
        #constitutive_driver
        t̂, ∂t∂Ĵ, _ = material_response(material, Ĵ, materialstate[qp])

        #if iszero(t̂)
        #    continue
        #end

        t = R ⋅ t̂
        ∂t∂J = R ⋅ ∂t∂Ĵ ⋅ R'

        for i in 1:ndofs
            J̇ = shape_value(cv, qp, i)
            fe[i] += (J ⋅ ∂t∂J ⋅ J̇) * dΓ #(t ⋅ J̇) * dΓ
        end
    end

end

function integrate_dissipation!(element::CohesiveElement{dim_s,CV}, 
    elementstate::Vector{<:AbstractElementState}, 
    material::AbstractMaterial, 
    materialstate::AbstractArray{<:AbstractMaterialState}, 
    fe::Vector{T}, 
    ge::Base.RefValue{T},
    coords, 
    Δue::Vector,
    ue::Vector,
    due::Vector,
    Δt::T) where {dim_s,CV,T}
    
    cv = element.cv
    ndofs = Ferrite.ndofs(element)
    xe = coords + reinterpret(Vec{dim_s,T}, ue)

    reinit!(cv, xe)

    for qp in 1:getnquadpoints(cv)
        
        #Rotation matrix
        R = getR(cv,qp)

        dΓ = getdetJdA(cv, qp) * element.thickness2d
        
        J = function_value(cv, qp, ue)
        Ĵ = R'⋅J

        #
        g, dgdĴ = constitutive_driver_dissipation(material, Ĵ, materialstate[qp])
        dgdJ =  R ⋅ dgdĴ
        #ddd = material isa MatCZBilinear ? 3 : 2
        #dgdJhat = Vec{2,T}((materialstate[qp].dgdJ[1], materialstate[qp].dgdJ[ddd]))
        #g = materialstate[qp].g
        #dgdJ =  R ⋅ dgdJhat


        ge[] += g * dΓ 
        for i in 1:ndofs
            δJ = shape_value(cv, qp, i)
            fe[i] += dgdJ ⋅ δJ * dΓ #(t ⋅ J̇) * dΓ
        end
    end

end

function integrate_forcevector!(element::CohesiveElement{dim_s}, 
    elementstate::Vector{<:AbstractElementState}, 
    material::AbstractMaterial, 
    materialstate::AbstractArray{<:AbstractMaterialState}, 
    stresses::Vector{<:SymmetricTensor{2,3,T}},
    strains::Vector{<:SymmetricTensor{2,3,T}},    
    fe::Vector{T}, 
    cell, 
    Δue::Vector,
    ue::Vector,
    due::Vector,
    Δt::T) where {dim_s,T}
    
    cv = element.cv
    
    ndofs = Ferrite.ndofs(element)

    xe = cell + reinterpret(Vec{dim_s,T}, ue)

    reinit!(cv, xe)
    
    for qp in 1:getnquadpoints(cv)
        
        #Rotation matrix
        R = getR(cv,qp)

        dΓ = getdetJdA(cv, qp) * element.thickness2d
        
        J = function_value(cv, qp, ue)
        Ĵ = R'⋅J
        
        #constitutive_driver
        t̂, _, new_matstate = material_response(material, Ĵ, materialstate[qp])
        materialstate[qp] = new_matstate

        t = R ⋅ t̂
        for i in 1:ndofs
            δui = shape_value(cv, qp, i)
            fe[i] += (t ⋅ δui) * dΓ
        end
    end

end

function integrate_massmatrix!(element::CohesiveElement, elementstate::Vector{<:AbstractElementState}, material::AbstractMaterial, cell, me::AbstractMatrix, ue::AbstractVector, due::AbstractVector)
    
end

function bodyforce!(element::CohesiveElement, elstate::AbstractElementState, material::AbstractMaterial, cell, fe::AbstractVector, forcevec::AbstractVector)
    
end