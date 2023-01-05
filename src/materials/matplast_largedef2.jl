export MatHyperElasticPlastic2, MatHyperElasticPlastic2State

"""
    MatHyperElasticPlastic2 <: AbstractMaterial

Hyper elastic plastic material for large deformations
"""
struct MatHyperElasticPlastic2{M <: HyperElasticMaterial} <: HyperElasticMaterial
    density::Float64
    elastic_material::M #Elastic material, Yeoh or Neo-hook
    τ₀::Float64		    #Yield stress
    H::Float64		    #Hardening
end

struct MatHyperElasticPlastic2State <: AbstractMaterialState
    S::SymmetricTensor{2,3,Float64,6} #2nd PK stress
    ϵᵖ::Float64					       #Plastic strain
    Fᵖ::Tensor{2,3,Float64,9}		   #Plastic deformation grad
    ν::Tensor{2,3,Float64,9}	
end

is_dissipative(::MatHyperElasticPlastic2) = true

# # # # # # #
# Constructors
# # # # # # #

function getmaterialstate(::MatHyperElasticPlastic2)
    S = zero(SymmetricTensor{2,3,Float64})
    Fᵖ = one(Tensor{2,3,Float64})
    ν = zero(Tensor{2,3,Float64})
    return MatHyperElasticPlastic2State(S, 0.0, Fᵖ, ν)
end

get_material_state_type(m::MatHyperElasticPlastic2) = MatHyperElasticPlastic2State

function MatHyperElasticPlastic2(; 
    elastic_material::HyperElasticMaterial,
    τ₀  ::T,
    H   ::T,
    density ::T = NaN, 
    ) where T

    return MatHyperElasticPlastic2(density, elastic_material, τ₀, H)
end


#Yield function Hyperelastic
function hyper_yield_function2(C::SymmetricTensor, Fᵖ::Tensor, ϵᵖ::T, emat::MatNeoHook, τ₀::Float64, H::Float64) where T
    Cᵉ = symmetric(transpose(Fᵖ)⋅C⋅Fᵖ)

    S̃, ∂S̃∂Cₑ = _constitutive_driver(emat, Cᵉ)

    #Mandel stress
    Mbar = Cᵉ⋅S̃
    Mᵈᵉᵛ = dev(Mbar)

    yield_func = sqrt(3/2)*norm(Mᵈᵉᵛ) - (τ₀ + H*ϵᵖ)

    return yield_func, Mᵈᵉᵛ, Cᵉ, S̃, ∂S̃∂Cₑ
end

function constitutive_driver(mp::MatHyperElasticPlastic2, C::SymmetricTensor{2,3}, state::MatHyperElasticPlastic2State)

    #Is it possible to combine the computation for S and ∂S∂E (using auto diff?)
    S, ∂S∂C, ϵᵖ, ν, Fᵖ, _, _ = _compute_2nd_PK_2(mp, C, state)
    
    #Numerical diff:
    # func(C) = _compute_2nd_PK(mp, C, state)[1]
    # ∂S∂C_num = numdiff(func, C)

    return S, ∂S∂C, MatHyperElasticPlastic2State(S, ϵᵖ, Fᵖ, ν)
end

function _compute_2nd_PK_2(mp::MatHyperElasticPlastic2, C::SymmetricTensor{2,dim,T}, state::MatHyperElasticPlastic2State) where {dim,T}
    
    Fᵖ = state.Fᵖ
    ϵᵖ = state.ϵᵖ
    Δγ = 0.0
    I = one(C)
    Iᵈᵉᵛ = otimesu(I,I) - 1/3*otimes(I,I)
    Cᵛ = tovoigt(Tensor{2,3}((i,j) -> C[i,j]))

    phi, Mᵈᵉᵛ, Cᵉ, S̃, ∂S̃∂Cₑ = hyper_yield_function2(C, state.Fᵖ, state.ϵᵖ, mp.elastic_material, mp.τ₀, mp.H)
    dFᵖdC = zero(Tensor{4,dim,T})
    dΔγdC = zero(Tensor{2,dim,T})
    dgdC = zero(Tensor{2,dim,T})
    g = 0.0

    if phi > 0
        residual(x) = compute_residual(x, Cᵛ, mp.H, mp.τ₀, mp.elastic_material, state.ϵᵖ, state.Fᵖ)
        jacoban(x) = ForwardDiff.jacobian(residual, x)
        x0 = vcat(tovoigt(state.Fᵖ), 0.0)

        res = nlsolve(residual, jacoban, x0; iterations=20, ftol = 1e-9, method=:newton, show_trace = false)
        x_res = res.zero

        RY(Cᵛ) = compute_residual(x_res, Cᵛ, mp.H, mp.τ₀, mp.elastic_material, state.ϵᵖ, state.Fᵖ)

        drdX = jacoban(x_res)
        drdY = ForwardDiff.jacobian(RY, Cᵛ)

        dXdY = -drdX \ drdY

        dFᵖdC = fromvoigt(Tensor{4,3}, dXdY[1:9,1:9])
        dΔγdC = fromvoigt(Tensor{2,3}, dXdY[10,:])

        Δγ = x_res[10]
        Fᵖ = fromvoigt(Tensor{2,3,Float64}, x_res[1:9])
        ϵᵖ += Δγ

        #Recompute
        Cᵉ = symmetric(transpose(Fᵖ)⋅C⋅Fᵖ)
        S̃, ∂S̃∂Cₑ = _constitutive_driver(mp.elastic_material, Cᵉ)

        Mbar = Cᵉ⋅S̃
        Mᵈᵉᵛ = dev(Mbar)
    end

    S = symmetric(Fᵖ ⋅ S̃ ⋅ transpose(Fᵖ))
    ν = sqrt(3/2)*Mᵈᵉᵛ/(norm(Mᵈᵉᵛ))

    ∂S∂Fᵖ = otimesu(I, (Fᵖ ⋅ S̃)) + otimesl((Fᵖ ⋅ S̃), I)
    ∂S∂S̃ = otimesu(Fᵖ,Fᵖ)
    dCₑdC = otimesu(transpose(Fᵖ),transpose(Fᵖ))
    ∂Cₑ∂Fᵖ = otimesl(I, transpose(Fᵖ) ⋅ C) + otimesu(transpose(Fᵖ) ⋅ C, I)
    dSdC = ∂S∂Fᵖ ⊡ dFᵖdC + ∂S∂S̃ ⊡ ∂S̃∂Cₑ ⊡ (dCₑdC + ∂Cₑ∂Fᵖ ⊡ dFᵖdC)

    if phi > 0.0
        ∂M∂Cₑ = otimesu(I, S̃)
        ∂M∂S̃ = otimesu(Cᵉ, I)
        M = Cᵉ ⋅ S̃
        g, dgdC = _compute_dissipation2(M, ν, Δγ, dCₑdC, ∂S̃∂Cₑ, dΔγdC, ∂M∂Cₑ, ∂M∂S̃, Iᵈᵉᵛ, I)
    end

    return S, dSdC, ϵᵖ, ν, Fᵖ, g, dgdC
end

function compute_residual(x, Cᵛ, H::Float64, τ₀::Float64, emat::HyperElasticMaterial, ⁿϵᵖ::Float64, ⁿFᵖ::Tensor)

    #Extract
    C = fromvoigt(Tensor{2,3}, Cᵛ)
    Fᵖ = fromvoigt(Tensor{2,3}, x[1:9])
    Δγ  = x[10]

    #Mandel stress
    Cᵉ = symmetric(transpose(Fᵖ)⋅C⋅Fᵖ)
    S̃, ∂S̃∂Cₑ = _constitutive_driver(emat, Cᵉ)
    Mbar = Cᵉ⋅S̃
    Mᵈᵉᵛ = dev(Mbar)
    ν = sqrt(3/2)*Mᵈᵉᵛ/(norm(Mᵈᵉᵛ))

    #Residuals
    R1 = Fᵖ - ⁿFᵖ + Δγ * ⁿFᵖ⋅ν
    R2 = √(3/2)*norm(Mᵈᵉᵛ) - (τ₀ + H*(ⁿϵᵖ + Δγ))

    return vcat(tovoigt(R1),R2)
end


function constitutive_driver_dissipation(mp::MatHyperElasticPlastic2, C::SymmetricTensor{2,3}, state::MatHyperElasticPlastic2State) 
    
    #Is it possible to combine the computation for S and ∂S∂E (using auto diff?)
    _, _, _, _, _, g, dgdC = _compute_2nd_PK_2(mp, C, state)

    #dgdC, g = gradient((C)->_compute_2nd_PK(mp, C, state)[6], C, :all)

    return g, dgdC
end


function _compute_dissipation2(M, ν, Δγ, dCₑdC, ∂S̃∂Cₑ, dΔγdC, ∂M∂Cₑ, ∂M∂S̃, Iᵈᵉᵛ, I)

    Mᵈᵉᵛ = dev(M)

    ∂g∂M = ν * Δγ
    ∂g∂ν = M * Δγ
    ∂g∂Δγ = M ⊡ ν
    
    ∂ν∂Mᵈᵉᵛ = √(3/2) * (1/norm(Mᵈᵉᵛ)) * (Iᵈᵉᵛ - (Mᵈᵉᵛ ⊗ Mᵈᵉᵛ)/(norm(Mᵈᵉᵛ)^2) )
    ∂Mᵈᵉᵛ∂M = Iᵈᵉᵛ
    dMdC = (∂M∂Cₑ + ∂M∂S̃ ⊡ ∂S̃∂Cₑ) ⊡ dCₑdC

    g = M ⊡ ν*Δγ
    dgdC = (∂g∂M + ∂g∂ν ⊡ ∂ν∂Mᵈᵉᵛ ⊡ ∂Mᵈᵉᵛ∂M) ⊡ dMdC  +  ∂g∂Δγ * dΔγdC 

    # ∂g∂Cₑ = otimesu(I, S̃) ⊡ (Δγ*ν)
    # ∂g∂S̃ = otimesu(Cₑ, I) ⊡ (Δγ*ν)
    #dgdC = ∂g∂Cₑ ⊡ dCₑdC  +  ∂g∂S̃ ⊡ ∂S̃∂Cₑ ⊡ dCₑdC   +   ∂g∂Δγ * dΔγdC   +  ∂g∂ν ⊡ ∂ν∂Mᵈᵉᵛ ⊡ ∂Mᵈᵉᵛ∂M ⊡ dMdC
    
    return g, dgdC
end