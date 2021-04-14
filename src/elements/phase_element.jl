export PhaseFieldElement
"""
PhaseFieldElement{dim,order,shape,T,CV<:JuAFEM.Values}

"""

struct PhaseFieldElement{dim,order,shape,T,CV1<:JuAFEM.Values,CV2<:JuAFEM.Values} <: AbstractElement
    thickness::T #used in 2d
    Gc::T
    lc::T
    μ::T
    λ::T

    celltype::Type{<:Cell}
    cv_u::CV1
    cv_d::CV2
    g::Function
end

struct PhaseFieldElementState <: AbstractElementState
    H::Vector{Float64}
end
PhaseFieldElementState(e::PhaseFieldElement, d::Float64 = 0.0) = PhaseFieldElementState(fill(d, getnquadpoints(e.cv_u)))
getelementstate(e::PhaseFieldElement) = PhaseFieldElementState(zeros(Float64, getnquadpoints(e.cv_u)))

get_elementstate_type(::PhaseFieldElement) = PhaseFieldElementState

JuAFEM.getnquadpoints(e::PhaseFieldElement) = getnquadpoints(e.cv_u)
JuAFEM.ndofs(e::PhaseFieldElement) = getnbasefunctions(e.cv_u) + getnbasefunctions(e.cv_d)
JuAFEM.getcelltype(e::PhaseFieldElement) = e.celltype
getncoords(s::PhaseFieldElement) = JuAFEM.getngeobasefunctions(s.cv_u)

has_constant_massmatrix(::PhaseFieldElement) = true

get_fields(::PhaseFieldElement{dim,order,shape,T}) where {dim,order,shape,T} = return [Field(:u, Lagrange{dim,shape,order}(), dim), Field(:d, Lagrange{dim,shape,order}(), 1)]

function PhaseFieldElement{dim,order,refshape,T}(;
    thickness = 1.0,     
    Gc::T,
    lc::T,
    μ::T,
    λ::T,
    qr_order::Int=2, 
    celltype::Type{<:Cell}) where {dim, order, refshape, T}
    
    qr = QuadratureRule{dim, refshape}(qr_order)
    
    ip = Lagrange{dim, refshape, order}()
    geom_ip = JuAFEM.default_interpolation(celltype)

    cv_u = CellVectorValues(qr, ip, geom_ip)
    cv_d = CellScalarValues(qr, ip, geom_ip)

    δ(i,j) = i == j ? 1.0 : 0.0
    g(i,j,k,l) = λ*δ(i,j)*δ(k,l) + μ*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))

    return PhaseFieldElement{dim,order,refshape,T,typeof(cv_u),typeof(cv_d)}(thickness, Gc, lc, μ, λ, celltype, cv_u, cv_d, g)
end

function integrate_forcevector_and_stiffnessmatrix!(element::PhaseFieldElement{dim, order, shape, T}, 
        elementstate::PhaseFieldElementState, 
        material::AbstractMaterial, 
        materialstate::AbstractVector{<:AbstractMaterialState},
        ke::Matrix,
        fe::Vector, 
        cell, 
        Δue::Vector,
        a::Vector,
        due::Vector,
        Δt::T) where {dim, order, shape, T}

        
    cv_u = element.cv_u
    cv_d = element.cv_d
    reinit!(cv_u, cell)
    reinit!(cv_d, cell)
    ndofs_u = getnbasefunctions(cv_u)
    ndofs_d = getnbasefunctions(cv_d)

    ue = a[1:ndofs_u]
    de = a[(1:ndofs_d) .+ ndofs_u]

    I = one(SymmetricTensor{2,dim,T})
    δ(i,j) = i == j ? 1.0 : 0.0
    I1 = SymmetricTensor{4,dim}((i,j,k,l)->δ(i,j)*δ(k,l))
    I2 = SymmetricTensor{4,dim}((i,j,k,l)->δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))

    λ = element.λ
    μ = element.μ
    lc = element.lc
    Gc = element.Gc

    for qp in 1:getnquadpoints(cv_u)
        
        ∇u = function_gradient(cv_u, qp, ue)
        ∇d = function_gradient(cv_d, qp, de)
        d = function_value(cv_d, qp, de)
        
        dΩ = getdetJdV(cv_u, qp) * element.thickness
       
        mac₊(x) = 0.5(x + abs(x))
        mac₋(x) = 0.5(x - abs(x))

        ε = symmetric(∇u)

        #p, n = eigen(ε)
        #ε⁺ = sum( i -> mac₊(p[i]) * (n[i,:] ⊗ n[i,:]), 1:dim) 
        #ε⁻ = sum( i -> mac₋(p[i]) * (n[i,:] ⊗ n[i,:]), 1:dim) 

        #σ⁺ = element.λ*tr(ε⁺)*I + 2*element.μ*ε⁺
        #σ⁻ = element.λ*tr(ε⁻)*I + 2*element.μ*ε⁻

        #σ = (1-d)^2 * σ⁺ + σ⁻ 
        #ψ⁺ = (1-d)^2 * (0.5*element.λ*mac₊(tr(ε))^2 + element.μ*tr(ε⁺ ⋅ ε⁺))
        #H = max(ψ⁺, elementstate.H[qp])
        #ΔH = ψ⁺ - H
        #∂H∂ε = ΔH >= 0.0 ? σ⁺ : zero(σ⁺)

        ψ    = (1-d)^2 * (0.5*λ*tr(ε)^2 + μ*ε⊡ε)
        σ    = (1-d)^2 * (λ*tr(ε)*I + 2*μ*ε)
        ∂σ∂ε = (1-d)^2 * (λ*I1 + μ*I2  )
        σ⁺ = σ

        H = max(ψ, elementstate.H[qp])
        ΔH = ψ - H
        ∂H∂ε = ΔH >= 0.0 ? σ : zero(σ)
        #∂H∂ε = σ
        elementstate.H[qp] = H

        # Hoist computations of δE
        for i in 1:ndofs_u
            ∇Ni = shape_symmetric_gradient(cv_u, qp, i)

            fe[i] += (σ ⊡ ∇Ni) * dΩ

            for j in 1:ndofs_u
                ∇Nj = shape_symmetric_gradient(cv_u, qp, j)

                ke[i,j] += (∇Ni ⊡ ∂σ∂ε ⊡ ∇Nj) *dΩ 
            end

            for j in 1:ndofs_d
                Nj = shape_value(cv_d, qp, j)
                
                ke[i,j+ndofs_u] += (2*(d-1)*Nj*(∇Ni ⊡ σ⁺)) *dΩ 
            end
        end

        for i in 1:ndofs_d
            Ni = shape_value(cv_d, qp, i)
            ∇Ni = shape_gradient(cv_d, qp, i)

            fe[i+ndofs_u] += (((Gc/lc + 2*H)*d - 2*H)*Ni + Gc*lc*∇d⋅∇Ni) * dΩ

            for j in 1:ndofs_u
                ∇Nj = shape_symmetric_gradient(cv_u, qp, j)

                ke[i + ndofs_u,j] += (2*(d-1)*Ni*(∇Nj ⊡ ∂H∂ε)) *dΩ 
            end

            for j in 1:ndofs_d
                Nj = shape_value(cv_d, qp, j)
                ∇Nj = shape_gradient(cv_d, qp, j)

                ke[i + ndofs_u, j + ndofs_u] += ((Gc/lc + 2*H)*Ni*Nj + Gc*lc*∇Ni⋅∇Nj) *dΩ 
            end

        end

    end

end


function integrate_dissipation!(
    element       :: PhaseFieldElement{dim, order, shape, T}, 
    elementstate  :: AbstractElementState,
    material      :: AbstractMaterial,
    materialstate :: AbstractArray{<:AbstractMaterialState},
    fe            :: Vector{T},
    ge            :: Base.RefValue{T},
    coords, 
    a             :: Vector,
    Δa            :: Vector,
    due           :: Vector,
    Δt            :: T
    ) where {dim, order, shape, T}
    
    cv_u = element.cv_u
    cv_d = element.cv_d
    reinit!(cv_d, coords)
    ndofs_u = getnbasefunctions(cv_u)
    ndofs_d = getnbasefunctions(cv_d)

    de = a[(1:ndofs_d) .+ ndofs_u]
    Δde = Δa[(1:ndofs_d) .+ ndofs_u]

    I = one(SymmetricTensor{2,dim,T})
    λ = element.λ
    μ = element.μ
    lc = element.lc
    Gc = element.Gc

    for qp in 1:getnquadpoints(cv_u)
        
        ∇d = function_gradient(cv_d, qp, de)
        d = function_value(cv_d, qp, de)

        Δd = function_value(cv_d, qp, Δde)
        Δ∇d = function_gradient(cv_d, qp, Δde)
        
        dΩ = getdetJdV(cv_u, qp) * element.thickness
       
        ge[] += (Gc/lc) * (d⋅Δd + lc^2*∇d⋅Δ∇d) * dΩ
        for i in 1:ndofs_d
            Ni = shape_value(cv_d, qp, i)
            ∇Ni = shape_gradient(cv_d, qp, i)

            fe[i+ndofs_u] += (Gc/lc) * (Ni*Δd + d*Ni + lc^2*(∇Ni⋅Δ∇d + ∇d⋅∇Ni)) * dΩ
        end

    end
end
