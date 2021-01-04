export MatVanDenBosch, MatVanDenBoschState

abstract type AbstractCohesiveMaterial <: AbstractMaterial end 

density(m::AbstractCohesiveMaterial) = 0.0

"""
interface_damage(m::AbstractCohesiveMaterial, dim::Int)

Returns the damage in the direction specifies by `dim`. 
"""
interface_damage(m::AbstractCohesiveMaterial, ::Int) = error("Not impelemented for cohesive material $m")

"""
max_traction_force(m::AbstractCohesiveMaterial, dim::Int)

Returns the maximum traction that can occur in the cohesive zone `dim`. 
ypicaly σₘₐₓ in the normal direction (3), and τₘₐₓ in the shearing directions (1 and 2).
"""
max_traction_force(m::AbstractCohesiveMaterial, ::Int) = error("Not impelemented for cohesive material $m")

"""
onset_displacement(m::AbstractCohesiveMaterial, dim::Int)

Returns the displacement/jump when the material starts to become damaged
"""
onset_displacement(m::AbstractCohesiveMaterial, ::Int) = error("Not impelemented for cohesive material $m")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Material X - Mixed Cohesive material
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


