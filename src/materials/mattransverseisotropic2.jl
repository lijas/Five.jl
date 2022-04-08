export MatTransvIsotropic, MatTransvIsotropicState

"""
    MatTransvIsotropic()
"""
struct MatTransvIsotropic <: MaterialModels.AbstractMaterial
    L⊥::Float64
    L₌::Float64
    M₌::Float64
    G⊥::Float64
    G₌::Float64
    density::Float64
end

function MatTransvIsotropicEngineeringConstants(; 
    E_L::T, 
    E_T::T, 
    G_LT::T,
    ν_LT::T, 
    ν_TT::T,
    ρ::T = 0.0) where T

    M₌ = (E_L^2*(ν_TT - 1))/(2*E_T*ν_LT^2 - E_L + E_L*ν_TT)
    L₌ = -(E_L*E_T*ν_LT)/(2*E_T*ν_LT^2 - E_L + E_L*ν_TT)
    L⊥ = -(E_T*(E_T*ν_LT^2 + E_L*ν_TT))/((ν_TT + 1)*(2*E_T*ν_LT^2 - E_L + E_L*ν_TT))
    G⊥ = E_T/(2*(ν_TT + 1))
    G₌ = G_LT

    return MatTransvIsotropic(L⊥, L₌, M₌, G⊥, G₌, ρ)
end

struct MatTransvIsotropicState <: MaterialModels.AbstractMaterialState
    σ::SymmetricTensor{2,3,Float64,6}
    a3::Vec{3,Float64}
end

# # # # # # #
# Constructors
# # # # # # #

function MaterialModels.initial_material_state(::MatTransvIsotropic, a3::Vec{3,Float64} = zero(Vec{3,Float64}))
    σ = zero(SymmetricTensor{2,3,Float64})
    return MatTransvIsotropicState(σ,a3)
end


# # # # # # #
# Drivers
# # # # # # #

function MaterialModels.material_response(m::MatTransvIsotropic, ε::SymmetricTensor{2,3,T,M}, state::MatTransvIsotropicState, Δt=nothing; cache=nothing, options=nothing) where {T,M}

    a3 = state.a3
    I = one(ε)
    A = symmetric( a3 ⊗ a3 )
    𝔸 = 0.25 * symmetric( otimesu(A,I) + otimesl(A,I) + otimesu(I,A) + otimesl(I,A) )
    Iˢʸᵐ = 0.5 * symmetric( (otimesu(I,I) + otimesl(I,I)) )
    
    #@assert( ismajorsymmetric(𝔸) )
    #@assert( isminorsymmetric(𝔸) )

    E = m.L⊥                                  * symmetric(I ⊗ I)     + 
        2m.G⊥                                 * Iˢʸᵐ                 + 
        (m.L₌ - m.L⊥)                         * symmetric(I⊗A + A⊗I) + 
        (m.M₌ - 4m.G₌ + 2m.G⊥ - 2m.L₌ + m.L⊥) * symmetric(A ⊗ A)     + 
        4(m.G₌ - m.G⊥)                        * 𝔸

    #@assert( ismajorsymmetric(E) )
    #@assert( isminorsymmetric(E) )
    σ = E ⊡ ε

    return σ, E, MatTransvIsotropicState(σ, state.a3)
end
