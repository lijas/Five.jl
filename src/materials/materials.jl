export LayeredMaterial, MatLinearElastic, Material2D
export PlaneStressMaterial, PlaneStrainMaterial

density(m::M) where M <: MaterialModels.AbstractMaterial = m.density

#@generated function materialstate_type(::Type{T}) where {T<:AbstractMaterial}
#
#end

materialstate_type(::Type{<:LinearElastic}) = LinearElasticState

"""
    is_dissipative
"""
is_dissipative(::AbstractMaterial) = false

function constitutive_driver_dissipation(::AbstractMaterial, ε::T, args...; kwargs...) where T
    return zero(Float64), zero(ε)
end

#include("matelastic.jl")
include("mattransvlinearelastic.jl")
include("mattransverseisotropic2.jl")
#include("matyeoh.jl")
#include("mathyper.jl")
#include("matplast_largedef.jl")


#include("cohesive/cohesive.jl")
include("cohesive/matczbilinear.jl")
#include("cohesive/matczbilinear_singlemode.jl")
include("cohesive/matczexponential.jl")
include("phasefield_material.jl")
#include("cohesive/mat_temp.jl")
#include("cohesive/matczbilinear2.jl")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# MatElasticSpring - Massless spring 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
"""
   MatElasticSpring
"""
struct MatElasticSpring <: MaterialModels.AbstractMaterial
    k::Float64
end

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Utility material for composites
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

struct LayeredMaterial{N,T,MAT} <: MaterialModels.AbstractMaterial
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
