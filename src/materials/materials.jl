"""
Abstract type for Materials
"""
abstract type AbstractMaterial end
abstract type HypoElasticMaterial <: AbstractMaterial end
abstract type HyperElasticMaterial <: AbstractMaterial end

abstract type AbstractMaterialState end

density(m::M) where M <: AbstractMaterial = m.density

"""
    constitutive_driver

The relation between the strains and the stress
"""
constitutive_driver

"""
    get_material_state_type

Get the MaterialState type for the material
"""
get_material_state_type


# # # #
# 2D materials
# # # #
@enum PLANE_STATE_2D plane_stress plane_strain

"""
Wrapper for making 3d-materials work in 3d
"""
struct Material2D{M <: AbstractMaterial}
    material::M
    plane_stress::PLANE_STATE_2D
end

include("matelastic.jl")
include("mattransvlinearelastic.jl")
include("matyeoh.jl")

include("cohesive/cohesive.jl")
include("cohesive/matczbilinear.jl")
include("cohesive/matczbilinear_singlemode.jl")
include("cohesive/matczexponential.jl")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# MatElasticSpring - Massless spring 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
"""
   MatElasticSpring
"""
struct MatElasticSpring <: AbstractMaterial
    k::Float64
end

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Utility material for composites
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

struct LayeredMaterial{N,T,MAT} <: AbstractMaterial
    materials::NTuple{N,MAT}
    orientations::NTuple{N,T}
end

function Base.getindex(m::LayeredMaterial, i::Int)
    return m.materials[i]
end

function LayeredMaterial(material::MAT, orientations::Vector{T}) where {T,MAT<: AbstractMaterial}
    nlayers = length(orientations)
    layermats = [deepcopy(material) for i in 1:nlayers]
    return LayeredMaterial{nlayers,T,MAT}(Tuple(layermats), Tuple(orientations))
end

function LayeredMaterial(layermats::Vector{MAT}, orientations::Vector{T}) where {T,MAT<: AbstractMaterial}
    @assert length(layermats) == length(orientations)
    nlayers = length(orientations)
    return LayeredMaterial{nlayers,T,MAT}(Tuple(layermats), Tuple(orientations))
end

nlayers(::LayeredMaterial{N}) where N = N

export LayeredMaterial, MatLinearElastic, Material2D