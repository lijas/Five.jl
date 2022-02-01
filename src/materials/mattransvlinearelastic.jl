export MatTransvLinearElastic, MatTransvLinearElasticState

"""
    MatTransvLinearElastic()

Standart material used for composites
"""
struct MatTransvLinearElastic <: MaterialModels.AbstractMaterial
    C::SymmetricTensor{4,3,Float64,36}
    α::Float64
    density::Float64
end

struct MatTransvLinearElasticState <: MaterialModels.AbstractMaterialState
    σ::SymmetricTensor{2,3,Float64,6}
end

# # # # # # #
# Constructors
# # # # # # #

function MaterialModels.initial_material_state(::MatTransvLinearElastic)
    σ = zero(SymmetricTensor{2,3,Float64})
    return MatTransvLinearElasticState(σ)
end

function MatTransvLinearElastic(;    
    E1::T,   
    E2::T,   
    E3::T = E2, 

    ν_12::T, 
    ν_23::T = ν_12, 
    ν_13::T = ν_12,  

    G_12::T, 
    G_23::T = G_12, 
    G_13::T = G_12,                           
    ρ::T = 0.0,
    α::T=0.0) where {T<:Real}

    #Complience matrix voight form
    C = [   1/E1 -ν_12/E1 -ν_13/E1   0       0         0;
        -ν_12/E1     1/E2 -ν_23/E2   0       0         0;
        -ν_13/E1 -ν_23/E2     1/E3   0       0         0;
            0         0        0     1/G_23  0         0;
            0         0        0     0       1/G_13    0;
            0         0        0     0       0       1/G_12];

    C = fromvoigt(SymmetricTensor{4,3}, inv(C))

    return MatTransvLinearElastic(C,α,ρ)
end


# # # # # # #
# Drivers
# # # # # # #
function _constitutive_driver(mp::MatTransvLinearElastic, ε::SymmetricTensor{2,3,T}) where T
    α = mp.α

    cα = cos(α)
    sα = sin(α)
    R = Tensor{2,3,T,9}((cα, sα, 0.0, -sα, cα, 0.0, 0.0, 0.0, 1.0))
    C = symmetric(otimesu(R,R) ⊡ mp.C ⊡ otimesu(R',R'))

    σ = C⊡ε

    return C, σ
end

function MaterialModels.material_response(mp::MatTransvLinearElastic, ε::SymmetricTensor{2,dim,T,M}, ::MatTransvLinearElasticState = getmaterialstate(mp)) where {dim,T,M}

    #dσ::SymmetricTensor{4,3,Float64,36}, σ::SymmetricTensor{2,3,Float64,6} = Tensors.gradient(e -> _constitutive_driver(mp, e), ε, :all)
    
    C, σ = _constitutive_driver(mp, ε)

    return σ, C, MatTransvLinearElasticState(σ)
end
