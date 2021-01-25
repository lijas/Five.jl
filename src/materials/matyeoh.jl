"""
    MaterialYeoh(; λ, μ, c_2, c_3)

Hyper elastic material, Yeoh
"""
struct MatYeoh <: HyperElasticMaterial
    density::Float64
    E::Float64
    nu::Float64 
    my::Float64
    lambda::Float64
end

struct MatYeohState <: AbstractMaterialState
    S::SymmetricTensor{2,3,Float64,6}
end

function getmaterialstate(::MatYeoh) 
    return MatYeohState()
end


get_material_state_type(m::MatYeoh) = MatYeohState
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

function getmaterialstate(::MatYeoh) 
    return MatYeohState()
end



function MatYeoh(; rho, E, nu)
    return MatYeoh(rho, E, nu, E / (2(1 + nu)), (E * nu) / ((1 + nu) * (1 - 2*nu)))
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
function constitutive_driver(mp::Union{MatYeoh, MatNeoHook}, C::SymmetricTensor{2,3}, state::AbstractMaterialState)
    ∂²Ψ∂C², ∂Ψ∂C, _Ψ =  hessian((C) -> ψ(mp, C), C, :all)
    S = 2.0 * ∂Ψ∂C
    ∂S∂C = 2.0 * ∂²Ψ∂C²
    return ∂S∂C, S, MaterialState(S)
end