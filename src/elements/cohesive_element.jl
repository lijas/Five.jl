export CohesiveElement

"""
CohesiveElement

Any order, any shape, any dim, 

"""

struct CohesiveElement{dim_p,dim_s,CV} <: AbstractElement
    thickness2d::Float64

    ndofs::Int
    fields::Vector{Field}
    celltype::Type{<:Cell}
    cv::CV
end

JuAFEM.getnquadpoints(e::CohesiveElement) = getnquadpoints(e.cv)
JuAFEM.ndofs(e::CohesiveElement) = e.ndofs
JuAFEM.getcelltype(e::CohesiveElement) = e.celltype

has_constant_massmatrix(::CohesiveElement) = true
getncoords(s::CohesiveElement) = JuAFEM.getngeobasefunctions(s.cv)

#TODO: Remove this constructor becuase it assumes CohesiveZone-type
function CohesiveElement{dim_s}(;order::Int,nqp::Int) where {dim_p,dim_s}

    ip = Lagrange{dim_s-1,RefCube,order}()
    mid_qr = QuadratureRule{dim_s-1,RefCube}(nqp)

    cv = SurfaceVectorValues(mid_qr,ip)

    nnodes = getnbasefunctions(ip)*2

    ndofs = nnodes*dim_s

    return CohesiveElement{dim_s-1,dim_s,typeof(cv)}(ndofs, [Field(:u, ip, dim_s)], Cell{dim_s,nnodes,1},cv)
end

function CohesiveElement{dim_p,dim_s,T}(cv::CV) where {dim_p,dim_s,T,CV}
    #@assert(dim_s == 2)
    order = 1

    nnodes = getnbasefunctions(cv) ÷ dim_s
    nfaces = 4

    ndofs = nnodes*dim_s

    ip = Lagrange{dim_s,RefCube,1}()

    return CohesiveElement{dim_p,dim_s,typeof(cv)}(ndofs, [Field(:u, ip, dim_s)], Cell{dim_s,nnodes,nfaces}, cv)    
end


function integrate_forcevector_and_stiffnessmatrix!(element::CohesiveElement{dim_p,dim_s,CV}, 
    elementstate::AbstractElementState, 
    material::AbstractMaterial, 
    materialstate::AbstractArray{<:AbstractMaterialState}, 
    ke::AbstractMatrix, 
    fe::Vector{T}, 
    cell, 
    Δue::Vector,
    ue::Vector,
    due::Vector,
    Δt::T) where {dim_p,dim_s,CV,T}
    
    cv = element.cv
    
    ndofs = JuAFEM.ndofs(element)

    xe = cell + reinterpret(Vec{dim_s,T}, ue)

    reinit!(cv, xe)

    A = 0
    for qp in 1:getnquadpoints(cv)
        
        #Rotation matrix
        R = getR(cv,qp)
        #@assert(det(R)==1.0)

        dΓ = getdetJdA(cv, qp) * info.thickness2d
        A += dΓ
        
        J = function_value(cv, qp, ue)
        Ĵ = R'⋅J
        
        #constitutive_driver
        t̂, ∂t∂Ĵ, new_matstate = constitutive_driver(material, Ĵ, materialstate[qp])
        materialstate[qp] = new_matstate

        #if iszero(t̂)
        #    continue
        #end

        t = R ⋅ t̂
        ∂t∂J = R ⋅ ∂t∂Ĵ ⋅ R'

        #
        if false
            ΔJ = function_value(cv, qp, Δue)
            K = 1.0
            ξ = 0.02*1
            σᵛ = ξ * K .* ΔJ/Δt
            ∂σᵛ∂J = ξ * K/Δt * one(SymmetricTensor{2,dim_s,T})

            t += σᵛ
            ∂t∂J += ∂σᵛ∂J
        end


        for i in 1:ndofs
            δui = shape_value(cv, qp, i)

            fe[i] += (t ⋅ δui) * dΓ
            for j in 1:ndofs
                δuj = shape_value(cv, qp, j)
                ke[i,j] += δui⋅∂t∂J⋅δuj * dΓ
            end
        end
    end

    return A
end

function integrate_fstar!(element::CohesiveElement{dim_p,dim_s,CV}, 
    elementstate::AbstractElementState, 
    material::AbstractMaterial, 
    materialstate::AbstractArray{<:AbstractMaterialState}, 
    fe::Vector{T}, 
    cell, 
    Δue::Vector,
    ue::Vector,
    due::Vector,
    Δt::T) where {dim_p,dim_s,CV,T}
    
    cv = element.cv
    
    ndofs = JuAFEM.ndofs(element)

    xe = cell + reinterpret(Vec{dim_s,T}, ue)

    reinit!(cv, xe)

    A = 0
    for qp in 1:getnquadpoints(cv)
        
        #Rotation matrix
        R = getR(cv,qp)
        #@assert(det(R)==1.0)

        dΓ = getdetJdA(cv, qp) * info.thickness2d
        A += dΓ
        
        J = function_value(cv, qp, ue)
        Ĵ = R'⋅J
        
        #constitutive_driver
        t̂, ∂t∂Ĵ, _ = constitutive_driver(material, Ĵ, materialstate[qp])

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

    return A
end

function integrate_forcevector!(element::CohesiveElement{dim_p,dim_s}, 
    elementstate::AbstractElementState, 
    material::AbstractMaterial, 
    materialstate::AbstractArray{<:AbstractMaterialState}, 
    new_materialstate::AbstractArray{<:AbstractMaterialState}, 
    fe::Vector{T}, 
    cell, 
    Δue::Vector,
    ue::Vector,
    due::Vector,
    Δt::T) where {dim_p,dim_s,T}
    
    cv = element.cv
    
    ndofs = JuAFEM.ndofs(element)

    xe = cell + reinterpret(Vec{dim_s,T}, ue)

    reinit!(cv, xe)
    
    A = 0
    for qp in 1:getnquadpoints(cv)
        
        #Rotation matrix
        R = getR(cv,qp)

        dΓ = getdetJdA(cv, qp) * info.thickness2d
        A += dΓ
        
        J = function_value(cv, qp, ue)
        Ĵ = R'⋅J
        
        #constitutive_driver
        t̂, _, new_matstate = constitutive_driver(material, Ĵ, materialstate[qp])
        new_materialstate[qp] = new_matstate

        if iszero(t̂)
            continue
        end
        t = R ⋅ t̂


        for i in 1:ndofs
            δui = shape_value(cv, qp, i)

            fe[i] += (t ⋅ δui) * dΓ

        end
    end

    return A
end

function integrate_massmatrix!(element::CohesiveElement, elementstate::AbstractElementState, material::AbstractMaterial, cell, me::AbstractMatrix, ue::AbstractVector, due::AbstractVector)
    
end

function bodyforce!(element::CohesiveElement{dim,T,M}, elstate::AbstractElementState, material::AbstractMaterial, cell, fe::AbstractVector, forcevec::AbstractVector) where {dim,T,M}
    
end