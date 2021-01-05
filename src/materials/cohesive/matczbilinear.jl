
export MatCZBilinear, MatCZBilinearState

"""
MatCZBilinear

Standard bi-linear cohesive zone law material model

From: 
Camanho, P., & Davila, C. G. (2002). Mixed-Mode Decohesion Finite Elements in for the Simulation Composite of Delamination Materials.

"""

struct MatCZBilinear{T} <: AbstractCohesiveMaterial
    K::T #Initial stiffness
    Gᴵ::NTuple{3,T} #Fracture toughness
    τᴹᵃˣ::NTuple{3,T} #Streanghts 
    δ⁰::NTuple{3,T} #Value of jump when traction is maximum
    δᶠ::NTuple{3,T} #Complete seperation
    η::T
end

struct MatCZBilinearState{T} <: AbstractMaterialState
    δᴹᵃˣₘ::T
    d::T
    t::Vec{3,T}
    J::Vec{3,T}
end

function MatCZBilinear(; K::T, Gᴵ::NTuple{3,T}, τᴹᵃˣ::NTuple{3,T}, η::T) where {T}
    δᶠ = 2 .*Gᴵ./τᴹᵃˣ
    δ⁰ = τᴹᵃˣ./K
    return MatCZBilinear{T}(K, Gᴵ, τᴹᵃˣ, δ⁰, δᶠ, η)
end

function getmaterialstate(mp::MatCZBilinear{T}, d::T=zero(T)) where {T}
    t = zero(Vec{3,T}) 
    J = zero(Vec{3,T})

    δᶠₘ = sqrt(mp.δᶠ[1]^2 + mp.δᶠ[2]^2)
    δ⁰ₘ = mp.δ⁰[1]
    δᴹᵃˣₘ = _delta_max(d, δᶠₘ, δ⁰ₘ)
    δᴹᵃˣₘ += δᴹᵃˣₘ*0.01
    return MatCZBilinearState(δᴹᵃˣₘ, d,t,J)
end

get_material_state_type(::MatCZBilinear{T}) where {T} = MatCZBilinearState{T}

function constitutive_driver(mp::MatCZBilinear{T}, J2d::Vec{2,T}, prev_state::MatCZBilinearState) where {T}
    #Pad with zero
    J = Vec{3,T}((J2d[1], zero(T), J2d[2]))

    #Call 3d routine
    t, δᴹᵃˣₘ, d = _constitutive_driver(mp, J, prev_state)
    dt::Tensor{2,3,T,9}, t::Vec{3,T} = JuAFEM.gradient(J -> _constitutive_driver(mp, J, prev_state)[1], J, :all)

    #Remove third direction
    t2d = Vec{2,T}((t[1], t[2]))
    dt2d = SymmetricTensor{2,2,T,3}((dt[1,1], dt[3,1], dt[3,3]))

    return t2d, dt2d, MatCZBilinearState(δᴹᵃˣₘ,d,t,J)
end

function constitutive_driver(mp::MatCZBilinear{T}, J::Vec{3,T}, prev_state::MatCZBilinearState) where {T}
    
    t, δᴹᵃˣₘ, d = _constitutive_driver(mp, J, prev_state)
    dt::Tensor{2,3,T,9}, t::Vec{3,T} = JuAFEM.gradient(J -> _constitutive_driver(mp, J, prev_state)[1], J, :all)
    
    return t, dt, MatCZBilinearState(δᴹᵃˣₘ,d,t,J)
end

#
_damage(δ, δᶠₘ, δ⁰ₘ) = (δ<=δ⁰ₘ) ? 0.0 : (δᶠₘ*(δ - δ⁰ₘ))/(δ*(δᶠₘ-δ⁰ₘ))
_delta_max(d, δᶠₘ, δ⁰ₘ) = (d==0.0) ? 0.0 : (-δᶠₘ*δ⁰ₘ)/(d*(δᶠₘ-δ⁰ₘ)-δᶠₘ)
interface_damage(state::MatCZBilinearState, ::Int = 1) = state.d
initial_stiffness(mat::MatCZBilinear) = mat.K

max_traction_force(mat::MatCZBilinear, dim::Int) = mat.τᴹᵃˣ[dim]
max_traction_force(mat::MatCZBilinear) = mat.τᴹᵃˣ

onset_displacement(mat::MatCZBilinear) = mat.δ⁰
onset_displacement(mat::MatCZBilinear, dim::Int) = mat.δ⁰[dim]

function _constitutive_driver(mp::MatCZBilinear{T1}, δ::Vec{dim,T2}, prev_state::MatCZBilinearState) where {dim,T1,T2}

    K = mp.K
    
    δᵢⱼ = one(Tensor{2,dim,T1})
    K == 0.0 && return zero(Vec{dim,T2}), 0.0, 1.0

    δˢʰᵉᵃʳ₀ = mp.δ⁰[1] # = δ⁰[2]

    macl(x) = (x<=0 ? (0.0) : x)

    #First check if it has failed
    if prev_state.d == 1.0
        #Only at traction force if negative penatration
        if δ[dim] <= 0.0
            D = δᵢⱼ[:,dim]⊗δᵢⱼ[dim,:] * K #* macl(-δ₃)/-δ₃
        else
            D = zero(Tensor{2,dim,T1})
        end
        return D⋅δ, prev_state.δᴹᵃˣₘ, prev_state.d
    end
    δˢʰᵉᵃʳ = dim==3 ? sqrt(δ[1]^2 + δ[2]^2) : δ[1]
    δₘ = sqrt((δˢʰᵉᵃʳ)^2 + macl(δ[dim])^2)

    δᴹᵃˣₘ = max(prev_state.δᴹᵃˣₘ, δₘ)

    β = δˢʰᵉᵃʳ/δ[dim]

    δ⁰ₘ = _onset_softening(δ[dim], mp.δ⁰, δˢʰᵉᵃʳ₀, β)
    δᶠₘ = _cohesive_bk_criterion(δ[dim], δ⁰ₘ, mp.δᶠ, K, β, mp.η, mp.Gᴵ)
    local D, d
    if δᴹᵃˣₘ <= δ⁰ₘ
        D = K*δᵢⱼ
        d=0.0
    elseif δ⁰ₘ < δᴹᵃˣₘ < δᶠₘ
        d = (δᶠₘ*(δᴹᵃˣₘ - δ⁰ₘ))/(δᴹᵃˣₘ*(δᶠₘ-δ⁰ₘ))
        D = δᵢⱼ * (1-d)*K
        if δ[dim] <= 0.0
            D += δᵢⱼ[:,dim]⊗δᵢⱼ[dim,:] * d * K# * macl(-δ₃)/-δ₃
        end
    elseif δᴹᵃˣₘ >= δᶠₘ
        D = zero(δᵢⱼ)
        if δ[dim] <= 0.0
            D += δᵢⱼ[:,dim]⊗δᵢⱼ[dim,:] * K# * macl(-δ₃)/-δ₃
        end
        d=1.0
    else
        error("MatCZBilinear failed")
    end
    
    return D⋅δ, δᴹᵃˣₘ, d
end


_cohesive_bk_criterion(δ₃, δ⁰ₘ, δᶠ::NTuple{3,T}, K::T, β, η::T, G::NTuple{3,T}) where {T} = _cohesive_bk_criterion(δ₃, δ⁰ₘ, δᶠ[1], δᶠ[2], K, β, η, G[1], G[2])
_cohesive_bk_criterion(δ₃, δ⁰ₘ, δᶠ::NTuple{2,T}, K::T, β, η::T, G::NTuple{2,T}) where {T} = _cohesive_bk_criterion(δ₃, δ⁰ₘ, δᶠ[1], 0.0,   K, β, η, G[1], G[2])
function _cohesive_bk_criterion(δ₃, δ⁰ₘ, δᶠ₁::T, δᶠ₂::T, K::T, β, η::T, Gᴵ::T, Gᴵᴵ::T) where T<:AbstractFloat
    δᶠₘ = 0.0
    if δ₃ > 0
        δᶠₘ = 2/(K*δ⁰ₘ)*(Gᴵ + (Gᴵᴵ - Gᴵ)*(β^2/(1+β^2))^η)
    else
        δᶠₘ = sqrt((δᶠ₁)^2 + (δᶠ₂)^2)
    end
    return δᶠₘ
end

function _cohesive_power_law()
    #TODO: implement from article
end

_onset_softening(δ₃, δ⁰::NTuple{3,T}, δˢʰᵉᵃʳ₀, β) where {T} = _onset_softening(δ₃, δ⁰[1], δ⁰[3], δˢʰᵉᵃʳ₀, β)
_onset_softening(δ₃, δ⁰::NTuple{2,T}, δˢʰᵉᵃʳ₀, β) where {T} = _onset_softening(δ₃, δ⁰[1], δ⁰[2], δˢʰᵉᵃʳ₀, β)
function _onset_softening(δ₃, δ⁰₁, δ⁰₃, δˢʰᵉᵃʳ₀, β)
    δ⁰ₘ = 0.0
    if 0.0 < δ₃
        δ⁰ₘ = δ⁰₃*δ⁰₁*sqrt((1 + β^2)/(δ⁰₁^2 + (β*δ⁰₃)^2))
    else
        δ⁰ₘ = δˢʰᵉᵃʳ₀ 
    end
    return δ⁰ₘ
end