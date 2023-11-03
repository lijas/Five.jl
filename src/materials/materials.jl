
density(m::M) where M <: MaterialModels.AbstractMaterial = m.density
mac₊(x::T) where T = x < 0.0 ? zero(T) : x 
mac₋(x::T) where T = x > 0.0 ? zero(T) : x 
heviside(x::T) where T = x > 0.0 ? one(T) : zero(T)

#@generated function materialstate_type(::Type{T}) where {T<:AbstractMaterial}
#
#end

get_material_state_type(::Type{<:PhaseFieldSpectralSplit}) = PhaseFieldSpectralSplitState
get_material_state_type(::Type{<:LinearElastic}) = LinearElasticState
get_material_state_type(::Type{<:TransverselyIsotropic}) = TransverselyIsotropicState
get_material_state_type(::Type{<:Plastic}) = PlasticState

"""
    is_dissipative
"""
is_dissipative(::AbstractMaterial) = false

function constitutive_driver_dissipation(::AbstractMaterial, ε::T, args...; kwargs...) where T
    return zero(Float64), zero(ε)
end

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
