export MatTransvLinearElastic, MatTransvLinearElasticState

"""
    MatTransvLinearElastic()

Standart material used for composites
"""
struct MatTransvLinearElastic <: AbstractMaterial
    C::SymmetricTensor{4,3,Float64,36}
    α::Float64
    density::Float64
end

struct MatTransvLinearElasticState <: AbstractMaterialState
    σ::SymmetricTensor{2,3,Float64,6}
end

# # # # # # #
# Constructors
# # # # # # #

function getmaterialstate(::MatTransvLinearElastic)
    σ = zero(SymmetricTensor{2,3,Float64})
    return MatTransvLinearElasticState(σ)
end

function MatTransvLinearElastic(;    E1::T,   E2::T,   E3::T = E2, 
                                        ν_12::T, ν_23::T = ν_12, ν_13::T = ν_12, 
                                        G_12::T, G_23::T = G_12, G_13::T = G_12, 
                                        ρ = 0.0,
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

get_material_state_type(::MatTransvLinearElastic) = MatTransvLinearElasticState

# # # # # # #
# Drivers
# # # # # # #
function _constitutive_driver(mp::MatTransvLinearElastic, ε) 
    α = mp.α

    _R = [cos(α) -sin(α) 0; 
          sin(α)  cos(α) 0; 
          0            0 1]
          
    R = Tensor{2,3}(Tuple(_R))
    C = (otimesu(R,R) ⊡ mp.C ⊡ otimesu(R',R'))

    σ = C⊡ε

    return symmetric(σ)
end

function constitutive_driver(mp::MatTransvLinearElastic, ε::SymmetricTensor{2,dim,T,M}, ::MatTransvLinearElasticState = getmaterialstate(mp)) where {dim,T,M}
    σ = _constitutive_driver(mp, ε) 
    dσ::SymmetricTensor{4,3,Float64,36}, σ::SymmetricTensor{2,3,Float64,6} = JuAFEM.gradient(e -> _constitutive_driver(mp, e), ε, :all)
    
    return σ, dσ, MatTransvLinearElasticState(σ)
end

function constitutive_driver_elastic(::MatTransvLinearElastic, ::SymmetricTensor{2,3}, ::MatTransvLinearElasticState = getmaterialstate(mp))

    return 0.0, zero(SymmetricTensor{2,3})

end