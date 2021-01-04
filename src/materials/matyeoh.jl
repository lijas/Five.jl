"""
    MatYeoh(...) <: HyperElasticMaterial

Hyperelastic material law, with the free energy:
    \\phi(C) = \\mu...
"""
struct MatYeoh <: HyperElasticMaterial
    density::Float64
    E::Float64
    nu::Float64 
    my::Float64
    lambda::Float64
end

struct MatYeohState <: AbstractMaterialState end

function MatYeohState(::MatYeoh) 
    return MatYeohState()
end

get_material_state_type(m::MatYeoh) = MatYeohState

function MatYeoh(; rho, E, nu)
    return MatYeoh(rho, E, nu, E / (2(1 + nu)), (E * nu) / ((1 + nu) * (1 - 2*nu)))
end

function _compute_2nd_PK(mp::MatYeoh, E)
    I = one(E)
    C = 2E + one(E)
    invC = inv(C)
    J = sqrt(det(C))
    return mp.my *(I - invC) + mp.lambda * log(J) * invC
end

function constitutive_driver(mp::MatYeoh, E, ::MatYeohState)
    ∂S∂E, SPK = gradient(E -> _compute_2nd_PK(mp, E), E, :all)
    return SPK, ∂S∂E, MatYeohState()
end
