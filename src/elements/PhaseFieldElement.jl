export PhaseFieldElement
"""
PhaseFieldElement{dim,order,shape,T,CV<:Ferrite.Values}

"""

struct PhaseFieldElement{dim,order,shape,T,CV1<:Ferrite.Values,CV2<:Ferrite.Values} <: AbstractElement
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
    H::Float64
end

initial_element_state(::PhaseFieldElement) = PhaseFieldElementState(0.0)

elementstate_type(::Type{<:PhaseFieldElement}) = PhaseFieldElementState

Ferrite.getnquadpoints(e::PhaseFieldElement) = getnquadpoints(e.cv_u)
Ferrite.ndofs(e::PhaseFieldElement) = getnbasefunctions(e.cv_u) + getnbasefunctions(e.cv_d)
Ferrite.getcelltype(e::PhaseFieldElement) = e.celltype
getncoords(s::PhaseFieldElement) = Ferrite.getngeobasefunctions(s.cv_u)

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
    geom_ip = Ferrite.default_interpolation(celltype)

    cv_u = CellVectorValues(qr, ip, geom_ip)
    cv_d = CellScalarValues(qr, ip, geom_ip)

    δ(i,j) = i == j ? 1.0 : 0.0
    g(i,j,k,l) = λ*δ(i,j)*δ(k,l) + μ*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))

    return PhaseFieldElement{dim,order,refshape,T,typeof(cv_u),typeof(cv_d)}(thickness, Gc, lc, μ, λ, celltype, cv_u, cv_d, g)
end

function _integrate_forcevector!(element::PhaseFieldElement{dim, order, shape, T2}, 
    elementstate::AbstractVector{PhaseFieldElementState}, 
    a::Vector{T}, cv_u, cv_d) where {dim, order, shape, T, T2}

    ndofs_u = getnbasefunctions(cv_u)
    ndofs_d = getnbasefunctions(cv_d)
    
    fe = zeros(T, ndofs_u+ndofs_d)
        
    ue = a[1:ndofs_u]
    de = a[(1:ndofs_d) .+ ndofs_u]

    I = one(SymmetricTensor{2,dim})
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

       # p, n = eigen(ε)
        #ε⁺ = sum( i -> mac₊(p[i]) * (n[i,:] ⊗ n[i,:]), 1:dim) 
       #ε⁻ = sum( i -> mac₋(p[i]) * (n[i,:] ⊗ n[i,:]), 1:dim) 

        #σ⁺ = element.λ*mac₊(tr(ε))*I + 2*element.μ*ε⁺
        #σ⁻ = element.λ*mac₋(tr(ε))*I + 2*element.μ*ε⁻

        #σ = (1-d)^2 * σ⁺ + σ⁻ 
        σ    = (1-d)^2 * (λ*tr(ε)*I + 2*μ*ε)
        #ψ⁺ = (0.5*element.λ*mac₊(tr(ε))^2 + element.μ*tr(ε⁺ ⋅ ε⁺))
        ψ⁺ = (0.5*element.λ*tr(ε)^2 + element.μ*tr(ε ⋅ ε))
        H = max(ψ⁺, elementstate[qp].H)

        # Hoist computations of δE
        for i in 1:ndofs_u
            ∇Ni = shape_symmetric_gradient(cv_u, qp, i)

            fe[i] += (σ ⊡ ∇Ni) * dΩ
        end

        for i in 1:ndofs_d
            Ni = shape_value(cv_d, qp, i)
            ∇Ni = shape_gradient(cv_d, qp, i)

            fe[i+ndofs_u] += (((Gc/lc + 2*H)*d - 2*H)*Ni + Gc*lc*∇d⋅∇Ni) * dΩ
        end

    end

    return fe
end

function integrate_forcevector_and_stiffnessmatrix2!(
        element::PhaseFieldElement{dim, order, shape, T}, 
        elementstate::AbstractVector{PhaseFieldElementState}, 
        material::AbstractMaterial, 
        materialstate::AbstractVector{<:AbstractMaterialState},
        stresses::Vector{<:SymmetricTensor{2,3,T}},
        strains::Vector{<:SymmetricTensor{2,3,T}},
        ke::Matrix,
        fe::Vector, 
        cell, 
        Δue::Vector,
        a::Vector,
        due::Vector,
        Δt::T ) where {dim, order, shape, T}

    cv_u = element.cv_u
    cv_d = element.cv_d
    reinit!(cv_u, cell)
    reinit!(cv_d, cell)

    residual(a) = _integrate_forcevector!(element, elementstate, a, cv_u, cv_d) 

    _fe = residual(a)
    _ke = ForwardDiff.jacobian(residual, a)

    fe .= _fe
    ke .= _ke

    ndofs_u = getnbasefunctions(cv_u)
    ndofs_d = getnbasefunctions(cv_d)
    ue = a[1:ndofs_u]
    de = a[(1:ndofs_d) .+ ndofs_u]
    #Must store state variable outside
    for qp in 1:getnquadpoints(cv_u)
        
        ∇u = function_gradient(cv_u, qp, ue)
        d = function_value(cv_d, qp, de)
        mac₊(x) = 0.5(x + abs(x))
        mac₋(x) = 0.5(x - abs(x))

        ε = symmetric(∇u)

        #p, n = eigen(ε)
        #ε⁺ = sum( i -> mac₊(p[i]) * (n[i,:] ⊗ n[i,:]), 1:dim) 

        #ψ⁺ = (0.5*element.λ*mac₊(tr(ε))^2 + element.μ*tr(ε⁺ ⋅ ε⁺))
        ψ⁺ = (0.5*element.λ*tr(ε)^2 + element.μ*tr(ε ⋅ ε))
        H = max(ψ⁺, elementstate[qp].H)
        elementstate[qp] = PhaseFieldElementState(H)
    end

end

function integrate_forcevector_and_stiffnessmatrix!(
        element::PhaseFieldElement{dim, order, shape, T}, 
        elementstate::AbstractVector{PhaseFieldElementState}, 
        material::AbstractMaterial, 
        materialstate::AbstractVector{<:AbstractMaterialState},
        stresses::Vector{<:SymmetricTensor{2,3,T}},
        strains::Vector{<:SymmetricTensor{2,3,T}},
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
        ε = symmetric(∇u)
        
        dΩ = getdetJdV(cv_u, qp) * element.thickness


        #σ, dσdε, σ⁺, Ψ⁺ = phasefield_response(ε, d, Gc, λ, μ)
        σ, dσdε, σ⁺, Ψ⁺ = phasefield_response_AD(ε, d, Gc, λ, μ)
        #σ, dσdεe, _ = material_response(PlaneStrain(), material, ε, initial_material_state(material))

        
        if iszero(ε)
            σ, dσdε, _ = material_response(PlaneStrain(), material, ε, initial_material_state(material))
            ∂H∂ε = zero(σ)
            σ⁺ = zero(σ)
            H = 0.0
        end

        #println("elastic")
        #display("text/plain", tovoigt(dσdεe))
        #println("britt")
        #display("text/plain", tovoigt(dσdε))
        #println("ad")
        #display("text/plain", tovoigt(dσdε))
        #@show Ψ⁺ > elementstate.H[qp]
        
        if Ψ⁺ > elementstate[qp].H
            ∂H∂ε = σ⁺
            H = Ψ⁺
            elementstate[qp] = PhaseFieldElementState(H)
        else
            ∂H∂ε = zero(σ⁺)
            H = elementstate[qp].H
        end

        for i in 1:ndofs_u
            ∇Ni = shape_symmetric_gradient(cv_u, qp, i)

            fe[i] += (σ ⊡ ∇Ni) * dΩ

            for j in 1:ndofs_u
                ∇Nj = shape_symmetric_gradient(cv_u, qp, j)

                ke[i,j] += (∇Ni ⊡ dσdε ⊡ ∇Nj) *dΩ 
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

                ke[i+ndofs_u,j] += (2*(d-1)*Ni*(∇Nj ⊡ ∂H∂ε)) *dΩ 
            end

            for j in 1:ndofs_d
                Nj = shape_value(cv_d, qp, j)
                ∇Nj = shape_gradient(cv_d, qp, j)

                ke[i+ndofs_u, j+ndofs_u] += ((Gc/lc + 2*H)*Ni*Nj + Gc*lc*∇Ni⋅∇Nj) *dΩ 
            end

        end

    end

end


function integrate_dissipation!(
    element       :: PhaseFieldElement{dim, order, shape, T}, 
    elementstate  ::AbstractVector{PhaseFieldElementState}, 
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

mac₊(x::T) where T = x < 0.0 ? 0.0 : x #0.5(x + abs(x))
mac₋(x::T) where T = x > 0.0 ? 0.0 : x #0.5(x - abs(x))
heaviside(x::T) where T = x > 0.0 ? one(T) : zero(T)

function phasefield_response_AD(ε::SymmetricTensor{2,dim}, φ, Gc, λ, μ) where {dim}
    ε_dual = Tensors._load(ε, nothing)
    _σ, _σ⁺, _ψ⁺ = _phasefield_response(ε_dual, φ, Gc, λ, μ)

    σ = Tensors._extract_value(_σ)
    dσdε = Tensors._extract_gradient(_σ, ε)
    σ⁺ = Tensors._extract_value(_σ⁺)
    Ψ⁺ = Tensors._extract_value(_ψ⁺)

    return σ, dσdε, σ⁺, Ψ⁺
end

function _phasefield_response2(ε::SymmetricTensor{2,dim,T}, φ, Gc, λ, μ) where {dim,T}
    I = one(SymmetricTensor{2,dim})
    σ₀  = (λ*tr(ε)*I + 2*μ*ε)
    σ  = (1-φ)^2 * σ₀
    ψ⁺ = (0.5*λ*tr(ε)^2 + μ*tr(ε ⋅ ε))
    return σ, σ₀, ψ⁺
end

function _phasefield_response(ε::SymmetricTensor{2,dim,T}, φ, Gc, λ, μ) where {dim,T}
    I = one(SymmetricTensor{2,dim})

    p, n = Tensors.eigen(ε)
    ε⁺ = zero(SymmetricTensor{2,dim,T})
    ε⁻ = zero(SymmetricTensor{2,dim,T})
    
    for i in 1:2 
        ε⁺ += mac₊(p[i]) * symmetric((n[:,i] ⊗ n[:,i]))
        ε⁻ += mac₋(p[i]) * symmetric((n[:,i] ⊗ n[:,i]))
    end

    σ⁺ = λ*mac₊(tr(ε))*I + 2*μ*ε⁺
    σ⁻ = λ*mac₋(tr(ε))*I + 2*μ*ε⁻

    σ = (1-φ)^2 * σ⁺ + σ⁻ 
    ψ⁺ = (0.5*λ*mac₊(tr(ε))^2 + μ*(ε⁺ ⊡ ε⁺))

    return σ, σ⁺, ψ⁺
end


function phasefield_response(ε::SymmetricTensor{2,dim}, φ, Gc, λ, μ) where {dim}

    I2 = one(SymmetricTensor{2,dim})
    I4 = one(SymmetricTensor{4,dim})
    I2oI2 = symmetric(I2⊗I2)

    gφ = (1-φ)^2

    I⁺ = positive_proj(ε)
    I⁻ = I4 - I⁺

    ε⁺ = I⁺ ⊡ ε
    ε⁻ = ε - ε⁺

    trε = tr(ε)

    σ⁺ = λ * mac₊(trε) * I2 + 2.0*Gc * ε⁺
    σ⁻ = λ * mac₋(trε) * I2 + 2.0*Gc * ε⁻
    σ = gφ*σ⁺ + σ⁻

    Ψ⁺ = 0.5*λ*mac₊(trε)^2 + Gc*(ε⁺⊡ε⁺)

    H1 = heaviside(trε)
    H2 = 1.0 - H1

    dσ = gφ * (2.0*Gc*I⁺ + λ*H1*(I2oI2)) + 
              (2.0*Gc*I⁻ + λ*H2*(I2oI2))

    return σ, dσ, σ⁺, Ψ⁺
end

function positive_proj(ε::SymmetricTensor{2,dim}) where dim

    I⁺ = zero(SymmetricTensor{4,dim})
    epos = zeros(dim)
    diag = zeros(dim)

    v, p = eigen(ε)
    for i in 1:dim
        epos[i] = mac₊(v[i])
        if v[i] > 0.0
            diag[i]=1.0
        end
    end

    for i in 1:dim
        Ma = p[:,i] ⊗ p[:,i]
        I⁺ += symmetric(Ma ⊗ Ma) * diag[i]
    end

    for a in 1:dim
        for b in 1:dim
            Ma = p[:,a] ⊗ p[:,a]
            Mb = p[:,b] ⊗ p[:,b]

            Gab = otimesu(Ma,Mb) + otimesl(Ma,Mb)
            Gba = otimesu(Mb,Ma) + otimesl(Mb,Ma)

            if abs(v[a]-v[b]) <= 1e-13
                θab = 0.5*(diag[a]+diag[b])/2.0
            else
                θab = 0.5*(epos[a]-epos[b])/(v[a]-v[b])
            end
            I⁺ += θab*symmetric(Gab+Gba)
        end
    end

    return I⁺
end
