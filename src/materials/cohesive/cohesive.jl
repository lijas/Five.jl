export MatCZKolluri, MatCZKolluriState
export max_traction_force, interface_damage, onset_displacement

abstract type AbstractCohesiveMaterial <: AbstractMaterial end 

density(m::AbstractCohesiveMaterial) = 0.0

"""
    interface_damage(m::AbstractCohesiveMaterial, dim::Int)

Returns the damage in the direction specifies by `dim`. 
"""
interface_damage

"""
    interface_damage(m::AbstractCohesiveMaterial, dim::Int)

Returns all damage parameters for the material. 
"""
n_damage_parameters

"""
    max_traction_force(m::AbstractCohesiveMaterial, dim::Int)

Returns the maximum traction that can occur in the cohesive zone `dim`. 
ypicaly σₘₐₓ in the normal direction (3), and τₘₐₓ in the shearing directions (1 and 2).
"""
max_traction_force

"""
    onset_displacement(m::AbstractCohesiveMaterial, dim::Int)

Returns the displacement/jump when the material starts to become damaged
"""
onset_displacement


