export MatHyperElasticPlastic, MatHyperElasticPlasticState

"""
    MatHyperElasticPlastic <: AbstractMaterial

Hyper elastic plastic material for large deformations
"""
struct MatHyperElasticPlastic{M <: HyperElasticMaterial} <: HyperElasticMaterial
    rho::Float64
    elastic_material::M #Elastic material, Yeoh or Neo-hook
    τ₀::Float64		    #Yield stress
    H::Float64		    #Hardening
end

struct MatHyperElasticPlasticState <: AbstractMaterialState
    S::SymmetricTensor{2,3,Float64,6} #2nd PK stress
    ϵᵖ::Float64					       #Plastic strain
    Fᵖ::Tensor{2,3,Float64,9}		   #Plastic deformation grad
    ν::Tensor{2,3,Float64,9}	
end

is_dissipative(::MatHyperElasticPlastic) = true

# # # # # # #
# Constructors
# # # # # # #

function getmaterialstate(::MatHyperElasticPlastic)
    S = zero(SymmetricTensor{2,3,Float64})
    Fᵖ = one(Tensor{2,3,Float64})
    ν = zero(Tensor{2,3,Float64})
    return MatHyperElasticPlasticState(S, 0.0, Fᵖ, ν)
end

get_material_state_type(m::MatHyperElasticPlastic) = MatHyperElasticPlasticState

function MatHyperElasticPlastic(; 
    elastic_material::HyperElasticMaterial,
    τ₀  ::T,
    H   ::T,
    rho ::T = NaN, 
    ) where T

    return MatHyperElasticPlastic(rho, elastic_material, τ₀, H)
end

# # # # # # #
# Drivers
# # # # # # #


# Constitutive driver for hyperelastic-plastic material
# Unicode are used for variable names, so the code should be self-explanitory
function _compute_2nd_PK(mp::MatHyperElasticPlastic, C::SymmetricTensor{2,dim,T}, state::MatHyperElasticPlasticState) where {dim,T}
    emat, τ₀, H = (mp.elastic_material, mp.τ₀, mp.H)
    
    I = one(C)
    Iᵈᵉᵛ = otimesu(I,I) - 1/3*otimes(I,I)
    
    ϵᵖ = state.ϵᵖ
    Fᵖ = state.Fᵖ
    ν = state.ν
    
    dγ = 0.0
    
    phi, Mᵈᵉᵛ, Cᵉ, S̃, ∂S̃∂Cₑ = hyper_yield_function(C, state.Fᵖ, state.ϵᵖ, emat, τ₀, H)
    dFᵖdC = zero(Tensor{4,dim,T})
    dΔγdC = zero(Tensor{2,dim,T})
    dgdC = zero(Tensor{2,dim,T})
    g = 0.0

    if phi > 0 #Plastic step
        #Setup newton varibles
        TOL = 1e-9
        newton_error, counter = TOL + 1, 0
        dγ = 0.0
        #ⁿν = sqrt(3/2)*Mᵈᵉᵛ/(norm(Mᵈᵉᵛ))
        local ∂ϕ∂Mᵈᵉᵛ, ∂Mᵈᵉᵛ∂M, ∂M∂Cₑ, dϕdΔγ, ∂ϕ∂M, ∂M∂S̃
        while true; counter +=1 

            Fᵖ = state.Fᵖ - dγ*state.Fᵖ⋅state.ν
            ϵᵖ = state.ϵᵖ + dγ

            R, Mᵈᵉᵛ, Cₑ, S̃, ∂S̃∂Cₑ = hyper_yield_function(C, Fᵖ, ϵᵖ, emat, τ₀, H) 
            
            #Numerical diff
            #h = 1e-6
            #J = (yield_function(C, state.Fᵖ - (dγ+h)*state.Fᵖ⋅state.ν, state.ϵᵖ + (dγ+h), μ, λ, τ₀, H)[1] - yield_function(C, state.Fᵖ - dγ*state.Fᵖ⋅state.ν, state.ϵᵖ + dγ, μ, λ, τ₀, H)[1]) / h

            ∂ϕ∂M = sqrt(3/2) * (Mᵈᵉᵛ/norm(Mᵈᵉᵛ)) #⊡ Iᵈᵉᵛ #  state.ν#
            ∂M∂Cₑ = otimesu(I, S̃)
            ∂Cₑ∂Fᵖ = otimesl(I, transpose(Fᵖ)⋅C) + otimesu(transpose(Fᵖ)⋅C, I)
            ∂Fᵖ∂Δγ = -state.Fᵖ⋅state.ν
            ∂M∂S̃ = otimesu(Cₑ, I)

            dϕdΔγ = ∂ϕ∂M ⊡ ((∂M∂Cₑ + ∂M∂S̃ ⊡ ∂S̃∂Cₑ) ⊡ ∂Cₑ∂Fᵖ ⊡ ∂Fᵖ∂Δγ) - H
             
            ddγ = dϕdΔγ\-R
            dγ += ddγ

            newton_error = norm(R)


            if newton_error < TOL
                break
            end

            if counter > 6
                #@warn("Could not find equilibrium in material")
                break
            end
        end
        
        ∂Cₑ∂C = otimesu(transpose(Fᵖ),transpose(Fᵖ))
        ∂ϕ∂C = ∂ϕ∂M ⊡ ((∂M∂Cₑ + ∂M∂S̃ ⊡ ∂S̃∂Cₑ) ⊡ ∂Cₑ∂C) 
        dΔγdC = -inv(dϕdΔγ) * ∂ϕ∂C
        dFᵖdC = -(state.Fᵖ⋅state.ν) ⊗ dΔγdC

    end

    
    ν = norm(Mᵈᵉᵛ)==0.0 ? zero(Mᵈᵉᵛ) : sqrt(3/2)*Mᵈᵉᵛ/(norm(Mᵈᵉᵛ))

    S = (Fᵖ ⋅ S̃ ⋅ transpose(Fᵖ))

    ∂S∂Fᵖ = otimesu(I, (Fᵖ ⋅ S̃)) + otimesl((Fᵖ ⋅ S̃), I)
    ∂S∂S̃ = otimesu(Fᵖ,Fᵖ)
    ∂Cₑ∂Fᵖ = otimesl(I, transpose(Fᵖ) ⋅ C) + otimesu(transpose(Fᵖ) ⋅ C, I)
    ∂Cₑ∂C = otimesu(transpose(Fᵖ),transpose(Fᵖ))
    dCₑdC = ∂Cₑ∂C + ∂Cₑ∂Fᵖ ⊡ dFᵖdC
    dSdC = ∂S∂Fᵖ ⊡ dFᵖdC + ∂S∂S̃ ⊡ ∂S̃∂Cₑ ⊡ dCₑdC

    if phi > 0.0
        M = Cᵉ ⋅ S̃
        g, dgdC = _compute_dissipation(M, ν, dγ, dCₑdC, ∂S̃∂Cₑ, dΔγdC, ∂M∂Cₑ, ∂M∂S̃, Iᵈᵉᵛ, I)
    end

    return symmetric(S), dSdC, ϵᵖ, ν, Fᵖ, g, dgdC
end

#Yield function Hyperelastic
function hyper_yield_function(C::SymmetricTensor, Fᵖ::Tensor, ϵᵖ::T, emat::MatNeoHook, τ₀::Float64, H::Float64) where T
    Cᵉ = symmetric(transpose(Fᵖ)⋅C⋅Fᵖ)

    S̃, ∂S̃∂Cₑ = _constitutive_driver(emat, Cᵉ)

    #Mandel stress
    Mbar = Cᵉ⋅S̃
    Mᵈᵉᵛ = dev(Mbar)

    yield_func = sqrt(3/2)*norm(Mᵈᵉᵛ) - (τ₀ + H*ϵᵖ)

    return yield_func, Mᵈᵉᵛ, Cᵉ, S̃, ∂S̃∂Cₑ
end

function constitutive_driver(mp::MatHyperElasticPlastic, C::SymmetricTensor{2,3}, state::MatHyperElasticPlasticState)

    #Is it possible to combine the computation for S and ∂S∂E (using auto diff?)
    S, ∂S∂C, ϵᵖ, ν, Fᵖ, _, _ = _compute_2nd_PK(mp, C, state)
    
    #Numerical diff:
    # func(C) = _compute_2nd_PK(mp, C, state)[1]
    # ∂S∂C_num = numdiff(func, C)

    return S, ∂S∂C, MatHyperElasticPlasticState(S, ϵᵖ, Fᵖ, ν)
end



#=
function _compute_2nd_PK_2(mp::MatHyperElasticPlastic, C::SymmetricTensor{2,dim,T}, state::MatHyperElasticPlasticState) where {dim,T}
    
    Fᵖ = state.Fᵖ
    ϵᵖ = state.ϵᵖ
    Δγ = 0.0
    I = one(C)
    Iᵈᵉᵛ = otimesu(I,I) - 1/3*otimes(I,I)
    Cᵛ = tovoigt(Tensor{2,3}((i,j) -> C[i,j]))

    phi, Mᵈᵉᵛ, Cᵉ, S̃, ∂S̃∂Cₑ = hyper_yield_function(C, state.Fᵖ, state.ϵᵖ, mp.elastic_material, mp.τ₀, mp.H)
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

        dFᵖdC = fromvoigt(Tensor{4,3,Float64}, dXdY)
        dΔγdC = fromvoigt(Tensor{2,3,Float64}, dXdY[10,:])

        Δγ = x_res[10]
        Fᵖ = fromvoigt(Tensor{2,3,Float64}, x_res)
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
        g, dgdC = _compute_dissipation(M, ν, dγ, dCₑdC, ∂S̃∂Cₑ, dΔγdC, ∂M∂Cₑ, ∂M∂S̃, Iᵈᵉᵛ, I)
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
    R1 = Fᵖ - ⁿFᵖ + Δγ* ⁿFᵖ⋅ν
    R2 = √(3/2)*norm(Mᵈᵉᵛ) - (τ₀ + H*(ⁿϵᵖ + Δγ))

    return vcat(tovoigt(R1),R2)
end
=#

function constitutive_driver_dissipation(mp::MatHyperElasticPlastic, C::SymmetricTensor{2,3}, state::MatHyperElasticPlasticState) 
    
    #Is it possible to combine the computation for S and ∂S∂E (using auto diff?)
    _, _, _, _, _, g, dgdC = _compute_2nd_PK(mp, C, state)

    #dgdC, g = gradient((C)->_compute_2nd_PK(mp, C, state)[6], C, :all)

    return g, dgdC
end


function _compute_dissipation(M, ν, Δγ, dCₑdC, ∂S̃∂Cₑ, dΔγdC, ∂M∂Cₑ, ∂M∂S̃, Iᵈᵉᵛ, I)

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