struct PhaseFieldSpectralSplit <: MaterialModels.AbstractMaterial
    Gc::Float64
    lc::Float64
    μ::Float64
    λ::Float64
end

struct PhaseFieldSpectralSplitState <: MaterialModels.AbstractMaterialState
    σ⁺::SymmetricTensor{2,3,Float64,6} 
    Ψ⁺::Float64
end

function MaterialModels.initial_material_state(::PhaseFieldSpectralSplit)
    return PhaseFieldSpectralSplitState(zero(SymmetricTensor{2,3}), 0.0)
end

function MaterialModels.initial_material_state(::PhaseFieldSpectralSplit, history::Float64)
    return PhaseFieldSpectralSplitState(zero(SymmetricTensor{2,3}), history)
end

function PhaseFieldSpectralSplit(; E::Float64, ν::Float64, Gc::Float64, lc::Float64) 
    μ = E / (2(1 + ν))
    λ = (E * ν) / ((1 + ν) * (1 - 2ν))
    return PhaseFieldSpectralSplit(Gc, lc, μ, λ)
end

materialstate_type(::Type{<:PhaseFieldSpectralSplit}) = PhaseFieldSpectralSplitState

function my_material_response(mp::PhaseFieldSpectralSplit, ε::SymmetricTensor{2,dim}, φ::AbstractFloat, state::PhaseFieldSpectralSplitState, Δt=nothing; cache=nothing, options=nothing) where dim
    (; Gc, λ, μ) = mp
    
    local σ, dσdε, σ⁺, Ψ⁺
    if iszero(ε) || (dim==3 && count(x->iszero(x), ε[[1,5,9]]) == 2) #If two of the diagonal terms are zero (or equal to each other), the gradient of the eigen values will be NaN
        δ(i,j) = i == j ? 1.0 : 0.0
        f = (i,j,k,l) -> λ*δ(i,j)*δ(k,l) + μ*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))    
        dσdε = SymmetricTensor{4, dim, Float64}(f)
        σ⁺ = σ = zero(ε)
        Ψ⁺ = 0.0
    else 
        ε_dual = Tensors._load(ε, nothing)
        _σ, _σ⁺, _ψ⁺ = _phasefield_response(ε_dual, φ, Gc, λ, μ)
        
        σ    = Tensors._extract_value(_σ)
        dσdε = Tensors._extract_gradient(_σ, ε)
        σ⁺   = Tensors._extract_value(_σ⁺)
        Ψ⁺   = Tensors._extract_value(_ψ⁺)
        
        #σ, dσdε, σ⁺, Ψ⁺ = phasefield_response2(ε, φ, Gc, λ, μ)
    end

    if dim == 2
        σ⁺ = MaterialModels.increase_dim(σ⁺)
    end

    return σ, dσdε, PhaseFieldSpectralSplitState(σ⁺, Ψ⁺)
end

#=function _phasefield_response2(ε::SymmetricTensor{2,dim,T}, φ, Gc, λ, μ) where {dim,T}
    I = one(SymmetricTensor{2,dim})
    σ₀  = (λ*tr(ε)*I + 2*μ*ε)
    σ  = (1-φ)^2 * σ₀
    ψ⁺ = (0.5*λ*tr(ε)^2 + μ*tr(ε ⋅ ε))
    return σ, σ₀, ψ⁺
end=#

function _phasefield_response(ε::SymmetricTensor{2,dim,T}, φ, Gc, λ, μ) where {dim,T}
    I = one(SymmetricTensor{2,dim})

    p, n = Tensors.eigen(ε)
    ε⁺ = zero(SymmetricTensor{2,dim,T})
    ε⁻ = zero(SymmetricTensor{2,dim,T})
    
    for i in 1:dim
        _n = n[:,i]
        non = symmetric(_n ⊗ _n)
        ε⁺ += mac₊(p[i]) * non
        ε⁻ += mac₋(p[i]) * non
    end
    
    trε = tr(ε)
    σ⁺ = λ*mac₊(trε)*I + 2μ*ε⁺
    σ⁻ = λ*mac₋(trε)*I + 2μ*ε⁻

    σ = (1-φ)^2 * σ⁺ + σ⁻ 
    ψ⁺ = 0.5*λ*mac₊(trε)^2 + μ*(ε⁺ ⊡ ε⁺)

    return σ, σ⁺, ψ⁺
end


function phasefield_response2(ε::SymmetricTensor{2,dim}, φ, Gc, λ, μ) where {dim}

    I2 = one(SymmetricTensor{2,dim})
    I4 = one(SymmetricTensor{4,dim})
    I2oI2 = symmetric(I2⊗I2)

    gφ = (1-φ)^2

    I⁺ = positive_proj(ε)
    I⁻ = I4 - I⁺

    ε⁺ = I⁺ ⊡ ε
    ε⁻ = ε - ε⁺

    trε = tr(ε)

    σ⁺ = λ * mac₊(trε) * I2 + 2μ * ε⁺
    σ⁻ = λ * mac₋(trε) * I2 + 2μ * ε⁻

    σ = gφ*σ⁺ + σ⁻

    Ψ⁺ = 0.5*λ*mac₊(trε)^2 + μ*(ε⁺ ⊡ ε⁺)

    H1 = heaviside(trε)
    H2 = 1.0 - H1

    dσ = gφ * (2μ*I⁺ + λ*H1*(I2oI2)) + 
              (2μ*I⁻ + λ*H2*(I2oI2))

    return σ, dσ, σ⁺, Ψ⁺
end

function positive_proj(ε::SymmetricTensor{2,dim}) where dim

    I⁺ = zero(SymmetricTensor{4,dim})

    v, p = eigen(ε)

    epos = Vec{dim}(i -> mac₊(v[i]))
    diag = Vec{dim}(i -> v[i] > 0.0 ? 1.0 : 0.0)

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


struct PhaseFieldNoSplit <: MaterialModels.AbstractMaterial
    Gc::Float64
    lc::Float64
    μ::Float64
    λ::Float64
end

struct PhaseFieldNoSplitState <: MaterialModels.AbstractMaterialState
    σ⁺::SymmetricTensor{2,3,Float64,6} 
    Ψ⁺::Float64
end

function MaterialModels.initial_material_state(::PhaseFieldNoSplit)
    return PhaseFieldNoSplitState(zero(SymmetricTensor{2,3}), 0.0)
end

function MaterialModels.initial_material_state(::PhaseFieldNoSplit, history::Float64)
    return PhaseFieldNoSplitState(zero(SymmetricTensor{2,3}), history)
end

function PhaseFieldNoSplit(; E::Float64, ν::Float64, Gc::Float64, lc::Float64) 
    μ = E / (2(1 + ν))
    λ = (E * ν) / ((1 + ν) * (1 - 2ν))
    return PhaseFieldNoSplit(Gc, lc, μ, λ)
end

materialstate_type(::Type{<:PhaseFieldNoSplit}) = PhaseFieldNoSplitState

function my_material_response(mp::PhaseFieldNoSplit, ε::SymmetricTensor{2,dim}, φ::AbstractFloat, state::PhaseFieldNoSplitState, Δt=nothing; cache=nothing, options=nothing) where dim
    (; Gc, λ, μ) = mp
    
    ε_dual = Tensors._load(ε, nothing)
    _σ, _σ⁺, _ψ⁺ = _phasefield_response_nosplit(ε_dual, φ, Gc, λ, μ)
    
    σ    = Tensors._extract_value(_σ)
    dσdε = Tensors._extract_gradient(_σ, ε)
    σ⁺   = Tensors._extract_value(_σ⁺)
    Ψ⁺   = Tensors._extract_value(_ψ⁺)

    if dim == 2
        σ⁺ = MaterialModels.increase_dim(σ⁺)
    end

    return σ, dσdε, PhaseFieldNoSplitState(σ⁺, Ψ⁺)
end

function _phasefield_response_nosplit(ε::SymmetricTensor{2,dim,T}, φ, Gc, λ, μ) where {dim,T}
    I = one(SymmetricTensor{2,dim})
    σ₀  = (λ*tr(ε)*I + 2*μ*ε)
    σ  = (1-φ)^2 * σ₀
    ψ⁺ = (0.5*λ*tr(ε)^2 + μ*tr(ε ⋅ ε))
    return σ, σ₀, ψ⁺
end