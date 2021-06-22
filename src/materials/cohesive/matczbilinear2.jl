
export MatCZBilinear2, MatCZBilinear2State

"""
MatCZBilinear2

Standard bi-linear cohesive zone law material model

From: 
Turon, A., Camanho, P. P., Costa, J., & Dávila, C. G. (2006). A damage model for the simulation of delamination in advanced composites under variable-mode loading. Mechanics of Materials, 38(11), 1072–1089. https://doi.org/10.1016/j.mechmat.2005.10.003

"""

struct MatCZBilinear2{T} <: AbstractCohesiveMaterial
    K::T #Initial stiffness
    σₘₐₓ::T
    τₘₐₓ::T
    
    #Onset damage
    δ⁰ₙ::T
    δ⁰ₜ::T

    #Failure
    δᶠₙ::T
    δᶠₜ::T

    Gᴵ::T
    Gᴵᴵ::T

    η::T
end

struct MatCZBilinear2State{T} <: AbstractMaterialState
    r::T
    d::T
    t::Vec{3,T}
    Δ::Vec{3,T}
    #Dissipation
    g::T
    dgdJ::Vec{3,T}
end

function MatCZBilinear2(; K::T, Gᴵ, Gᴵᴵ, σₘₐₓ, τₘₐₓ, η::T) where {T}
    δᶠₙ = 2 *Gᴵ/σₘₐₓ
    δᶠₜ = 2 *Gᴵᴵ/τₘₐₓ
    δ⁰ₙ = σₘₐₓ/K
    δ⁰ₜ = τₘₐₓ/K
    return MatCZBilinear2{T}(K, σₘₐₓ, τₘₐₓ, δ⁰ₙ, δ⁰ₜ, δᶠₙ, δᶠₜ, Gᴵ, Gᴵᴵ, η)
end

function getmaterialstate(mp::MatCZBilinear2{T}, d::T=zero(T)) where {T}
    λ = _delta_max2(d, mp.δᶠₙ, mp.δ⁰ₙ)
    return MatCZBilinear2State(λ, d, zero(Vec{3,T}), zero(Vec{3,T}), 0.0, zero(Vec{3,T}))
end

get_material_state_type(::MatCZBilinear2{T}) where {T} = MatCZBilinear2State{T}

is_dissipative(::MatCZBilinear2) = true

function constitutive_driver(mp::MatCZBilinear2{T}, Δ::Vec{3,T}, prev_state::MatCZBilinear2State) where {T}
    
    t, r, d, g = _constitutive_driver(mp, Δ, prev_state)
    #dt::Tensor{2,3,T,9}, t::Vec{3,T} = Tensors.gradient(Δ -> _constitutive_driver(mp, Δ, prev_state)[1], Δ, :all)
    
    #dgdJ::Vec{3,T}, g = Tensors.gradient(Δ -> _constitutive_driver(mp, Δ, prev_state)[4], Δ, :all)
    
    dgdJ = zero(Vec{3,T})
    g = 0.0
    dt = zero(Tensor{2,3,T})
    return t, dt, MatCZBilinear2State(r, d, t, Δ, g, dgdJ)
end

function constitutive_driver_dissipation(mp::MatCZBilinear2{T}, J::Vec{3,T}, prev_state::MatCZBilinear2State) where {T}
    
    J_dual = Tensors._load(J, nothing)
    _, _, _, g_res = _constitutive_driver(mp, J_dual, prev_state)

    return Tensors._extract_value(g_res), Tensors._extract_gradient(g_res, J)
end

#
_damage2(δ, δᶠₘ, δ⁰ₘ) = (δ<=δ⁰ₘ) ? 0.0 : (δᶠₘ*(δ - δ⁰ₘ))/(δ*(δᶠₘ-δ⁰ₘ))
_delta_max2(d, δᶠₘ, δ⁰ₘ) = (d==0.0) ? δ⁰ₘ : (-δᶠₘ*δ⁰ₘ)/(d*(δᶠₘ-δ⁰ₘ)-δᶠₘ)
n_damage_parameters(state::MatCZBilinear2) = 1
interface_damage(state::MatCZBilinear2State, ::Int = 1) = state.d
initial_stiffness(mat::MatCZBilinear2) = mat.K

max_traction_force(mat::MatCZBilinear2, dim::Int) = mat.τᴹᵃˣ[dim]
max_traction_force(mat::MatCZBilinear2) = mat.τᴹᵃˣ

onset_displacement(mat::MatCZBilinear2) = mat.δ⁰
onset_displacement(mat::MatCZBilinear2, dim::Int) = mat.δ⁰[dim]

function _constitutive_driver(mp::MatCZBilinear2{T1}, Δ::Vec{3,T2}, state::MatCZBilinear2State) where {T1,T2}

    K = mp.K
    Δ⁰ₜ = mp.δ⁰ₜ
    δᵢⱼ = one(Tensor{2,3,T1})
    δ₃₃ = δᵢⱼ[:,3] ⊗ δᵢⱼ[:,3]
    d = state.d
    D0 = δᵢⱼ*K
    ΔΔ = Δ - state.Δ

    macl(x) = (x<=0 ? (0.0) : x)

    Δˢʰᵉᵃʳ = sqrt(Δ[1]^2 + Δ[2]^2)
    λ  = sqrt((Δˢʰᵉᵃʳ)^2 + macl(Δ[3])^2)

    r = max(state.r, λ)

    β = Δˢʰᵉᵃʳ/(Δˢʰᵉᵃʳ+macl(Δ[3]))
    isnan(β) && return K*δᵢⱼ⋅Δ, r, d, 0.0 

    B = β^2/(1 + 2β^2 - 2β)

    Δ⁰ = _onset_softening2(mp.δ⁰ₙ, mp.δ⁰ₜ, B, mp.η)
    Δᶠ = _cohesive_bk_criterion2(Δ⁰, mp.δ⁰ₙ, mp.δᶠₙ, mp.δ⁰ₜ, mp.δᶠₜ, B, mp.η)


    d = _damage2(r, Δᶠ, Δ⁰)

    t = (1-d)*D0⋅Δ - d*K*δᵢⱼ[:,3]*macl(-Δ[3]) 
    g =  0.5(Δ ⋅ (K*δᵢⱼ) ⋅ Δ) * d

    return t, r, d, g
end

function _cohesive_bk_criterion2( Δ⁰, Δ⁰₃, Δᶠ₃, Δ⁰ₜ, Δᶠₜ, B, η)
    return (Δ⁰₃*Δᶠ₃ + (Δ⁰ₜ*Δᶠₜ - Δ⁰₃*Δᶠ₃) * B^η)/Δ⁰
end

function _onset_softening2(Δ⁰ₙ, Δ⁰ₜ, B, η)
    return sqrt(Δ⁰ₙ^2 + ((Δ⁰ₜ)^2 - (Δ⁰ₙ)^2)*B^η)
end