export MatYeoh, MatNeoHook

"""
    MaterialYeoh(; λ, μ, c_2, c_3)

Hyper elastic material, Yeoh
"""
struct MatYeoh <: HyperElasticMaterial
    λ::Float64
    μ::Float64
    c_2::Float64
    c_3::Float64
end

struct MatYeohState <: AbstractMaterialState
    S::SymmetricTensor{2,3,Float64,6}
end

get_material_state_type(m::MatYeoh) = MatYeohState

function getmaterialstate(::MatYeoh)
    S = zero(SymmetricTensor{2,3,Float64})
    return MatYeohState(S)
end

function MatYeoh(; E::T, ν::T, c_2::T, c_3::T) where T <: AbstractFloat
    λ = (E*nu) / ((1+nu) * (1 - 2nu))
    μ = E / (2(1+nu))

    return MatYeoh(λ, μ, c_2, c_3)
end

"""
    MaterialNeoHook(; λ, μ)

Hyper elastic material, Neo-Hook
"""
struct MatNeoHook <: HyperElasticMaterial
    λ::Float64
    μ::Float64
end

struct MatNeoHookState <: AbstractMaterialState
    S::SymmetricTensor{2,3,Float64,6}
end

get_material_state_type(m::MatNeoHook) = MatNeoHookState

function getmaterialstate(::MatNeoHook)
    S = zero(SymmetricTensor{2,3,Float64})
    return MatNeoHookState(S)
end

function MatNeoHook(; E::T, ν::T) where T <: AbstractFloat
    λ = (E*ν) / ((1+ν) * (1 - ν))
    μ = E / (2(1+ν))

    return MatNeoHook(λ, μ)
end

#Free energy
function ψ(mp::MatYeoh, C::SymmetricTensor{2,3})
    J = sqrt(det(C))
    I = tr(C)
    return mp.μ/2 * (I-3) + mp.c_2*(I-3)^2 + mp.c_3*(I-3)^3 - mp.μ*log(J) + mp.λ/2 * log(J)^2
end

function ψ(mp::MatNeoHook, C::SymmetricTensor{2,3})
    J = sqrt(det(C))
    I = tr(C)
    return mp.μ/2 * (I-3) - mp.μ*log(J) + mp.λ/2 * log(J)^2
end

#Calculates the 2nd PK and C using Autodiff.
function constitutive_driver(mp::Union{MatYeoh, MatNeoHook}, C::SymmetricTensor{2,3}, ::AbstractMaterialState = getmaterialstate(mp))
    ∂²Ψ∂C², ∂Ψ∂C, _Ψ =  hessian((C) -> ψ(mp, C), C, :all)
    S = 2.0 * ∂Ψ∂C
    ∂S∂C = 2.0 * ∂²Ψ∂C²
    state = mp isa MatYeoh ? MatYeohState(S) : MatNeoHookState(S)
    return S, ∂S∂C, state
end