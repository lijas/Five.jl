"""
    MatLinearElastic() <: AbstractMaterial

Standar elastic material routine. 
"""
struct MatLinearElastic <: AbstractMaterial
    denisty::Float64
    E::Float64
    nu::Float64
    μ::Float64
    λ::Float64
    plain_stress::Bool
end

struct MatLinearElasticState <: AbstractMaterialState
    σ::SymmetricTensor{2,3,Float64,6}
end

# # # # # # #
# Constructors
# # # # # # #

MatLinearElasticState() = MatLinearElasticState(zero(SymmetricTensor{2,3,Float64}))

function MatLinearElasticState(::MatLinearElastic)
    σ = zero(SymmetricTensor{2,3,Float64})
    return MatLinearElasticState(σ)
end

get_material_state_type(m::MatLinearElastic) = MatLinearElasticState

function MatLinearElastic(; rho::T = NaN, E::T, nu::T, plane_stress::Bool=true) where T
    λ = (E*nu) / ((1+nu) * (1 - 2nu))
    μ = E / (2(1+nu))

    return MatLinearElastic(rho, E, nu, μ, λ, plane_stress)
end

# # # # # # #
# Drivers
# # # # # # #

function _constitutive_driver(mp::MatLinearElastic, ε::SymmetricTensor{2,3}) 
    σ = 2*mp.μ*ε + mp.λ*tr(ε)*one(SymmetricTensor{2,3,Float64})
    return σ
end

function _constitutive_driver(mp::MatLinearElastic, ε::SymmetricTensor{2,2}) 

    s = mp.plain_stress == true ? 1+mp.nu : 1/(1-mp.nu)
    σ = 2*mp.μ*(ε + ((s-1)/(2-s)) * tr(ε)*one(SymmetricTensor{2,2,Float64}))

    return σ
end

function constitutive_driver(mp::MatLinearElastic, ε::SymmetricTensor{2,3}, ::MatLinearElasticState = MatLinearElasticState())
    dσ,σ = JuAFEM.gradient(e -> _constitutive_driver(mp, e), ε, :all)
    return σ, dσ, MatLinearElasticState{3,Float64}(σ)
end

function constitutive_driver(mp::MatLinearElastic, ε::SymmetricTensor{2,2}, ::MatLinearElasticState = MatLinearElasticState()) 
    dσ, σ = JuAFEM.gradient(e -> _constitutive_driver(mp, e), ε, :all)

    #The 2d material should store state tensor in 3d
    σ33 = mp.plain_stress == true ? 0.0 : (mp.λ/(2*(mp.μ + mp.λ))) * (σ[1,1] + σ[2,2])
    σ_state = SymmetricTensor{2,3,Float64,6}((σ[1,1], 0.0, σ[1,2], σ33, 0.0, σ[2,2]))

    return σ, dσ, MatLinearElasticState(σ_state)
end

# 1d for bars
function constitutive_driver(mp::MatLinearElastic, ε::SymmetricTensor{2,1}, ::MatLinearElasticState = MatLinearElasticState()) 
    σ = mp.E*ε[1,1]

    #The 2d material should store state tensor in 3d
    σ_state = SymmetricTensor{2,3,Float64,6}((σ[1,1], 0.0, 0.0, 0.0, 0.0, 0.0))

    return σ, mp.E, MatLinearElasticState(σ_state)
end