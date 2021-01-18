export MatMixedCohesive, MatMixedCohesiveState
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Material X - Mixed Cohesive material
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

struct MatMixedCohesive{dim,T} <: AbstractCohesiveMaterial
    K::T #Initial stiffness
    Gᴵ::NTuple{dim,T} #Fracture toughness
    τᴹᵃˣ::NTuple{dim,T} #Streanghts 
    δ⁰::NTuple{dim,T} #Value of jump when traction is maximum
    δᶠ::NTuple{dim,T} #Complete seperation
    η::T
end

struct MatMixedCohesiveState{dim,T} <: AbstractMaterialState
    δᴹᵃˣₘ::T
    d::T
    t::Vec{dim,T}
    J::Vec{dim,T}

    #Dissipation
    g::T
    dgdJ::Vec{dim,T}
end

function MatMixedCohesive(K::T, Gᴵ::NTuple{dim,T}, τᴹᵃˣ::NTuple{dim,T}, η::T) where {dim,T}
    δᶠ = 2 .*Gᴵ./τᴹᵃˣ
    δ⁰ = τᴹᵃˣ./K
    return MatMixedCohesive{dim,T}(K, Gᴵ, τᴹᵃˣ, δ⁰, δᶠ, η)
end

function getmaterialstate(mp::MatMixedCohesive{dim,T}, d::T=zero(T)) where {dim,T}
    t = zero(Vec{dim,T}) #Traction vector
    J = zero(Vec{dim,T})
    #_cohesive_bk_criterion(δ₃, δ⁰ₘ, δᶠ₁, δᶠ₂, K, β, η, Gᴵ, Gᴵᴵ)
    δᶠₘ = dim==3 ? sqrt(mp.δᶠ[1]^2 + mp.δᶠ[2]^2) : mp.δᶠ[1]
    δ⁰ₘ = mp.δ⁰[1]
    δᴹᵃˣₘ = _delta_max(d, δᶠₘ, δ⁰ₘ)
    δᴹᵃˣₘ += δᴹᵃˣₘ*0.01
    t = zero(Vec{dim,T})
    J = zero(Vec{dim,T})
    return MatMixedCohesiveState(δᴹᵃˣₘ, d,t,J, zero(T), zero(Vec{dim,T}))
end

get_material_state_type(::MatMixedCohesive{dim,T}) where {dim,T} = MatMixedCohesiveState{dim,T}

interface_damage(state::MatMixedCohesiveState) = state.d

initial_stiffness(mat::MatMixedCohesive) = mat.K

max_traction_force(mat::MatMixedCohesive, dim::Int) = mat.τᴹᵃˣ[dim]
max_traction_force(mat::MatMixedCohesive) = mat.τᴹᵃˣ

onset_displacement(mat::MatMixedCohesive) = mat.δ⁰
onset_displacement(mat::MatMixedCohesive, dim::Int) = mat.δ⁰[dim]


function _constitutive_driver(mp::MatMixedCohesive{dim,T1}, δ::Vec{dim,T2}, prev_state::MatMixedCohesiveState) where {dim,T1,T2}

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
        return D⋅δ, prev_state.δᴹᵃˣₘ, prev_state.d, 0.0
    end
    δˢʰᵉᵃʳ = dim==3 ? sqrt(δ[1]^2 + δ[2]^2) : δ[1]
    δₘ = sqrt((δˢʰᵉᵃʳ)^2 + macl(δ[dim])^2)

    δᴹᵃˣₘ = max(prev_state.δᴹᵃˣₘ, δₘ)

    β = δˢʰᵉᵃʳ/δ[dim]

    δ⁰ₘ = _onset_softening(δ[dim], mp.δ⁰, δˢʰᵉᵃʳ₀, β)
    δᶠₘ = _cohesive_bk_criterion(δ[dim], δ⁰ₘ, mp.δᶠ, K, β, mp.η, mp.Gᴵ)
    local D, d, g
    if δᴹᵃˣₘ <= δ⁰ₘ
        D = K*δᵢⱼ
        d=0.0
        g=0.0
    elseif δ⁰ₘ < δᴹᵃˣₘ < δᶠₘ
        d = (δᶠₘ*(δᴹᵃˣₘ - δ⁰ₘ))/(δᴹᵃˣₘ*(δᶠₘ-δ⁰ₘ))
        D = δᵢⱼ * (1-d)*K
        Dg = K*δᵢⱼ # For compuation of dissipation
        if δ[dim] <= 0.0
            D += δᵢⱼ[:,dim]⊗δᵢⱼ[dim,:] * d * K# * macl(-δ₃)/-δ₃
            Dg -= δᵢⱼ[:,dim]⊗δᵢⱼ[dim,:] * K
        end

        Δd = d - prev_state.d
        g =  0.5(δ ⋅ Dg ⋅ δ) * Δd

    elseif δᴹᵃˣₘ >= δᶠₘ
        D = zero(δᵢⱼ)
        Dg = K*δᵢⱼ
        if δ[dim] <= 0.0
            D += δᵢⱼ[:,dim]⊗δᵢⱼ[dim,:] * K# * macl(-δ₃)/-δ₃
            Dg -= δᵢⱼ[:,dim]⊗δᵢⱼ[dim,:] * K
        end
        d=1.0
        Δd = d - prev_state.d
        g =  0.5(δ ⋅ Dg ⋅ δ) * Δd
    else
        @show δᴹᵃˣₘ, δ⁰ₘ, δ[dim], mp.δ⁰, δˢʰᵉᵃʳ₀, β, δᶠₘ
        error("Wrong D")
    end

    
    return D⋅δ, δᴹᵃˣₘ, d, g
end

function constitutive_driver(mp::MatMixedCohesive{dim,T}, J, prev_state::MatMixedCohesiveState) where {dim,T}
    
    t, δᴹᵃˣₘ, d, g = _constitutive_driver(mp, J, prev_state)
    #dt::Tensor{2,dim,T,M2},     t::Vec{dim,T} = JuAFEM.gradient(J -> _constitutive_driver(mp, J, prev_state)[1], J, :all)
    dt::Tensor{2,dim,T,dim^2}, t::Vec{dim,T} = JuAFEM.gradient(J -> _constitutive_driver(mp, J, prev_state)[1], J, :all)
    
    dgdJ::Vec{dim,T}, g = JuAFEM.gradient(J -> _constitutive_driver(mp, J, prev_state)[4], J, :all)
    
    return t, dt, MatMixedCohesiveState(δᴹᵃˣₘ,d,t,J, g, dgdJ)
end