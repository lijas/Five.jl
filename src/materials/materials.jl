export LayeredMaterial, MatLinearElastic, Material2D
export PlaneStressMaterial, PlaneStrainMaterial

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

"""
    getmaterialstate(::AbstractMaterial)

Construct a materialstate based in the input material
"""
getmaterialstate

"""
    is_dissipative
"""
is_dissipative(::AbstractMaterial) = false

include("matelastic.jl")
include("mattransvlinearelastic.jl")
#include("matyeoh.jl")
include("mathyper.jl")
include("matplast_largedef.jl")


include("cohesive/cohesive.jl")
include("cohesive/matczbilinear.jl")
include("cohesive/matczbilinear_singlemode.jl")
include("cohesive/matczexponential.jl")
#include("cohesive/mat_temp.jl")
include("cohesive/matczbilinear2.jl")

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

# # # #
# 2D materials
# # # #

@enum PLANE_STATE_2D PLANE_STRESS PLANE_STRAIN

"""
Wrapper for making 3d-materials work in 2d
"""
struct Material2D{M <: AbstractMaterial} <: AbstractMaterial
    material::M
    plane_state::PLANE_STATE_2D
end

PlaneStressMaterial(mat) = Material2D(mat, PLANE_STRESS)
PlaneStrainMaterial(mat) = Material2D(mat, PLANE_STRAIN)

is_dissipative(m::Material2D) = is_dissipative(m.material)

function getmaterialstate(m::Material2D)
    return getmaterialstate(m.material)
end

function constitutive_driver(m::Material2D{<:HyperElasticMaterial}, C::SymmetricTensor{2,2,T}, state::AbstractMaterialState) where T
    @assert(m.plane_state == PLANE_STRAIN)

    #Convert to 3d
    C₃ = SymmetricTensor{2,3,T,6}((C[1,1], zero(T), C[1,2], one(T), zero(T), C[2,2]))
    S₃, ∂S∂C₃, newstate₃ = constitutive_driver(m.material, C₃, state)

    #Convert back to 3d
    S = SymmetricTensor{2,2,T,3}((S₃[1,1],S₃[1,3],S₃[3,3]))
    ∂S∂C = SymmetricTensor{4,2,T,9}((∂S∂C₃[1,1,1,1], ∂S∂C₃[3,1,1,1], ∂S∂C₃[3,3,1,1], ∂S∂C₃[1,1,3,1], ∂S∂C₃[3,1,3,1], ∂S∂C₃[3,3,3,1], ∂S∂C₃[1,1,3,3], ∂S∂C₃[3,1,3,3], ∂S∂C₃[3,3,3,3]))

    return S, ∂S∂C, newstate₃
end

function constitutive_driver_dissipation(m::Material2D{<:HyperElasticMaterial}, C::SymmetricTensor{2,2,T}, state::AbstractMaterialState) where T
    @assert(m.plane_state == PLANE_STRAIN)

    #Convert to 3d
    C₃ = SymmetricTensor{2,3,T,6}((C[1,1], zero(T), C[1,2], one(T), zero(T), C[2,2]))
    g, dgdC₃ = constitutive_driver_dissipation(m.material, C₃, state)

    #Convert back to 3d
    dgdC = SymmetricTensor{2,2,T,3}((dgdC₃[1,1],dgdC₃[1,3],dgdC₃[3,3]))

    return g, dgdC 
end

function constitutive_driver(m::Material2D, ε_2d::SymmetricTensor{2,2,T}, prev_state::AbstractMaterialState) where T

    ε = SymmetricTensor{2,3,T,6}((ε_2d[1,1], T(0.0), ε_2d[1,2], T(0.0), T(0.0), ε_2d[2,2]))
    #ε = SymmetricTensor{2,3,T,6}((ε_2d[1,1], ε_2d[1,2], T(0.0), ε_2d[2,2], T(0.0), T(0.0)))

    local σ::SymmetricTensor{2,3,Float64,6}
    local C::SymmetricTensor{4,3,Float64,36}
    local state
    if m.plane_state == PLANE_STRESS
        #Newton variables
        NEWTON_TOL = 1e-8
        newton_counter = 0

        #Newton loop
        index = [2,4,6]
        εⱽ = tovoigt(ε)
        _error = NEWTON_TOL + 1.0
        σⱽ = zeros(Float64, 6)
        Cⱽ = zeros(Float64, 6, 6)
        while _error > NEWTON_TOL; 
            newton_counter +=1
            _σⱽ::SymmetricTensor{2,3,Float64,6}, _Cⱽ::SymmetricTensor{4,3,Float64,36}, state = constitutive_driver(m.material, ε, prev_state)
            tovoigt!(σⱽ, _σⱽ)
            tovoigt!(Cⱽ, _Cⱽ)

            Δεⱽ = -Cⱽ[index, index]\σⱽ[index]
            εⱽ[index] += Δεⱽ
            ε = fromvoigt(SymmetricTensor{2,3,T}, εⱽ)

            _error = norm(σⱽ[index])
            newton_counter == 5 && error("Solution not found in material")
        end

        #Modify stiffness matrix
        index2 = [1,3,5]
        Cmod = zeros(3,3)
        for (i,I) in enumerate(index2)
            for (j,J) in enumerate(index2)
                for B in index
                    Cmod[i,j] += (Cⱽ[I, B] * Cⱽ[B, J]) ./ Cⱽ[B,B]
                end
            end
        end
        Cⱽ[index2,index2] -= Cmod

        σ = fromvoigt(SymmetricTensor{2,3,T}, σⱽ)
        C = fromvoigt(SymmetricTensor{4,3,T}, Cⱽ)
    else
        σ, C, state = constitutive_driver(m.material, ε, prev_state)
    end

    σ_2d = SymmetricTensor{2,2,T,3}((σ[1,1],σ[1,3],σ[3,3]))
    C_2d = SymmetricTensor{4,2,T,9}((C[1,1,1,1], C[3,1,1,1], C[3,3,1,1], C[1,1,3,1], C[3,1,3,1], C[3,3,3,1], C[1,1,3,3], C[3,1,3,3], C[3,3,3,3]))

    #σ_2d = SymmetricTensor{2,2,T,3}((σ[1,1],σ[1,2],σ[2,2]))
    #C_2d = SymmetricTensor{4,2,T,9}((C[1,1,1,1], C[2,1,1,1], C[2,2,1,1], C[1,1,2,1], C[2,1,2,1], C[2,2,2,1], C[1,1,2,2], C[2,1,2,2], C[2,2,2,2]))
    
    return σ_2d, C_2d, state
end

get_material_state_type(m::Material2D) = get_material_state_type(m.material)

density(m::Material2D) = m.material.density