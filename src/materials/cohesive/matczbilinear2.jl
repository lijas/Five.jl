
export MatCZBilinear2, MatCZBilinear2State

"""
MatCZBilinear2

Standard bi-linear cohesive zone law material model

From: 
Camanho, P., & Davila, C. G. (2002). Mixed-Mode Decohesion Finite Elements in for the Simulation Composite of Delamination Materials.

"""

struct MatCZBilinear2{T} <: AbstractCohesiveMaterial
    K::T #Initial stiffness
    Gᴵ::NTuple{3,T} #Fracture toughness
    τᴹᵃˣ::NTuple{3,T} #Streanghts 
    δ⁰::NTuple{3,T} #Value of jump when traction is maximum
    δᶠ::NTuple{3,T} #Complete seperation
    η::T
end

struct MatCZBilinear2State{T} <: AbstractMaterialState
    δᴹᵃˣₘ::T
    d::T
    t::Vec{2,T}
    J::Vec{2,T}

    #Dissipation
    g::T
    dgdJ::Vec{2,T}
end

function MatCZBilinear2(; K::T, Gᴵ::NTuple{3,T}, τᴹᵃˣ::NTuple{3,T}, η::T) where {T}
    δᶠ = 2 .*Gᴵ./τᴹᵃˣ
    δ⁰ = τᴹᵃˣ./K
    return MatCZBilinear2{T}(K, Gᴵ, τᴹᵃˣ, δ⁰, δᶠ, η)
end

function getmaterialstate(mp::MatCZBilinear2{T}, d::T=zero(T)) where {T}
    t = zero(Vec{2,T}) 
    J = zero(Vec{2,T})

    δᶠₘ = mp.δᶠ[1]
    δ⁰ₘ = mp.δ⁰[1]
    δᴹᵃˣₘ = _delta_max(d, δᶠₘ, δ⁰ₘ)
    δᴹᵃˣₘ += δᴹᵃˣₘ*0.01
    return MatCZBilinear2State(δᴹᵃˣₘ, d,t,J, 0.0, zero(Vec{2,T}))
end

get_material_state_type(::MatCZBilinear2{T}) where {T} = MatCZBilinear2State{T}

#=function constitutive_driver(mp::MatCZBilinear2{T}, J::Vec{3,T}, prev_state::MatCZBilinear2State) where {T}
    
    t, δᴹᵃˣₘ, d, _ = _constitutive_driver(mp, J, prev_state)
    dt::Tensor{2,3,T,9}, t::Vec{3,T} = Tensors.gradient(J -> _constitutive_driver(mp, J, prev_state)[1], J, :all)
    
    dgdJ::Vec{dim,T}, g = Tensors.gradient(J -> _constitutive_driver(mp, J, prev_state)[4], J, :all)
    

    return t, dt, MatCZBilinear2State(δᴹᵃˣₘ,d,t,J, g, dgdJ)
end

function constitutive_driver_dissipation(mp::MatCZBilinear2{T}, J::Vec{3,T}, prev_state::MatCZBilinear2State) where {T}
    
    J_dual = Tensors._load(J)
    _, _, _, g_res = _constitutive_driver(mp, J_dual, prev_state)

    return Tensors._extract_value(g_res), Tensors._extract_gradient(g_res, J)
end

#2D
function constitutive_driver(mp::MatCZBilinear2{T}, J2d::Vec{2,T}, prev_state::MatCZBilinear2State) where {T}
    
    #Pad with zero
    J = Vec{3,T}((J2d[1], zero(T), J2d[2]))

    #Call 3d routine
    t, δᴹᵃˣₘ, d, _ = _constitutive_driver(mp, J, prev_state)
    dt::Tensor{2,3,T,9}, t::Vec{3,T} = Tensors.gradient(J -> _constitutive_driver(mp, J, prev_state)[1], J, :all)

    #Remove third direction
    t2d = Vec{2,T}((t[1], t[3]))
    dt2d = SymmetricTensor{2,2,T,3}((dt[1,1], dt[3,1], dt[3,3]))

    dgdJ::Vec{3,T}, g = Tensors.gradient(J -> _constitutive_driver(mp, J, prev_state)[4], J, :all)

    return t2d, dt2d, MatCZBilinear2State(δᴹᵃˣₘ,d,t,J, g, dgdJ)
end

function constitutive_driver_dissipation(mp::MatCZBilinear2{T}, J2d::Vec{2,T}, prev_state::MatCZBilinear2State) where {T}
    
    #Pad with zero
    J = Vec{3,T}((J2d[1], zero(T), J2d[2]))

    J_dual = Tensors._load(J)
    _, _, _, g_res = _constitutive_driver(mp, J_dual, prev_state)

    g, dg =  Tensors._extract_value(g_res), Tensors._extract_gradient(g_res, J)

    return g, Vec{2,T}((dg[1], dg[3]))
end=#


#=function constitutive_driver_dissipation(mp::MatCZBilinear2{T}, J, prev_state::MatCZBilinear2State) where {T}
    dim = 2
    t, δᴹᵃˣₘ, d, g = _constitutive_driver(mp, J, prev_state)
    #dt::Tensor{2,dim,T,M2},     t::Vec{dim,T} = Tensors.gradient(J -> _constitutive_driver(mp, J, prev_state)[1], J, :all)
    dt::Tensor{2,dim,T,dim^2}, t::Vec{dim,T} = Tensors.gradient(J -> _constitutive_driver(mp, J, prev_state)[1], J, :all)
    
    dgdJ::Vec{dim,T}, g = Tensors.gradient(J -> _constitutive_driver(mp, J, prev_state)[4], J, :all)
    
    return g, dgdJ
end=#


function constitutive_driver(mp::MatCZBilinear2{T}, J, prev_state) where {T}
    
    t, δᴹᵃˣₘ, d, g = _constitutive_driver(mp, J, prev_state)
    #dt::Tensor{2,dim,T,M2},     t::Vec{dim,T} = Tensors.gradient(J -> _constitutive_driver(mp, J, prev_state)[1], J, :all)
    dt::Tensor{2,2,T,2^2}, t::Vec{2,T} = Tensors.gradient(J -> _constitutive_driver(mp, J, prev_state)[1], J, :all)
    
    dgdJ::Vec{2,T}, g = Tensors.gradient(J -> _constitutive_driver(mp, J, prev_state)[4], J, :all)

    return t, dt, MatCZBilinear2State(δᴹᵃˣₘ,d, t, J, g, dgdJ)
end

#
interface_damage(state::MatCZBilinear2State, ::Int = 1) = state.d
initial_stiffness(mat::MatCZBilinear2) = mat.K

max_traction_force(mat::MatCZBilinear2, dim::Int) = mat.τᴹᵃˣ[dim]
max_traction_force(mat::MatCZBilinear2) = mat.τᴹᵃˣ

onset_displacement(mat::MatCZBilinear2) = mat.δ⁰
onset_displacement(mat::MatCZBilinear2, dim::Int) = mat.δ⁰[dim]

function _constitutive_driver(mp::MatCZBilinear2{T1}, δ::Vec{dim,T2}, prev_state::MatCZBilinear2State) where {dim,T1,T2}

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
