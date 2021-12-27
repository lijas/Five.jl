module Five

using Reexport
using ForwardDiff
using NLsolve

using Parameters
using TimerOutputs
using StaticArrays

using LinearAlgebra
using SparseArrays
using Statistics: mean
import Base: isdone

import Logging
import ProgressMeter
using InteractiveUtils

#import JLD2
import FileIO: save
using JLD2
#import IGA #https://github.com/lijas/IGA.jl.git


@reexport using Ferrite
@reexport using Tensors

LOG = Logging.global_logger()

abstract type AbstractPart{dim} end
abstract type AbstractPartState end

@enum SolverMode MODE1 MODE2
@enum SimulationTermination NORMAL_TERMINATION ERROR_TERMINATION ABORTED_SIMULATION _NO_TERMINATION

abstract type AbstractSolver{T} end

export SystemArrays, StateVariables, GlobalData

abstract type AbstractSystemArrays{T} end

"""
    SystemArrays(T::Type, ndofs::Int)

Contains global arrays, such as internal force vector and stiffness matrix.

**Possible changes**
Differnet solvers need different global arrays. For example, an static solver does need mass matrices,
while a explicit solver does. Maybe create different <: SystemArrays depending on solver. 

"""
mutable struct SystemArrays{T} <: AbstractSystemArrays{T}
    fⁱ::Vector{T}
    Kⁱ::SparseArrays.SparseMatrixCSC{T,Int}

    fᵉ::Vector{T} #External force vector
    Kᵉ::SparseArrays.SparseMatrixCSC{T,Int}
    
    Mᵈⁱᵃᵍ::SparseArrays.SparseMatrixCSC{T,Int}
    M    ::SparseArrays.SparseMatrixCSC{T,Int}
    
    #For arc-solvers
    q::Vector{T}

    #For dissipation solver:
    fᴬ::Vector{T}

    G::Base.RefValue{T}
end

function SystemArrays(T::Type, ndofs::Int)
    Mᵈⁱᵃᵍ = spzeros(T,ndofs,ndofs)
    M = spzeros(T,ndofs,ndofs)
    return SystemArrays(zeros(T,ndofs), spzeros(T,ndofs,ndofs), zeros(T,ndofs), spzeros(T,ndofs,ndofs), Mᵈⁱᵃᵍ, M, zeros(T,ndofs), zeros(T,ndofs), Ref(0.0))
end

abstract type AbstractStateVariables{T} end

"""
    StateVariables(T::Type, ndofs::Int)

Contains all information about the state of system (displacements, material damage etc).

**Field variables:**
* `d`: displacements 
* `v`: velocityes
* `a`: accelerations 
* `t`: current time
* `λ`: current loading factor (fᵉ = λ*f̂, used in arc-length solver) 
* `L`: Current step length (used in arc-length solvers)
* `Δd`, `Δv`, `Δa`, `Δt`, `Δλ`, `ΔL` - Difference between current and previous timestep of variables above

* `system_arrays`: Instance of [`SystemArrays`](@ref)

* `partstates`: Contains [`PartState`](@ref), one for each cell
* `prev_partstates`: Contains [`PartState`](@ref) from previous timesteps, one for each cell

* `step`: Number of steps taken up until this point

**Possible changes**
Remove all Δ-variables, and require two states instead, (previous and current) 

"""
mutable struct StateVariables{T} <: AbstractStateVariables{T}
    
    d::Vector{T}
    v::Vector{T}
    a::Vector{T}
    t::T
    Δt::T
    λ::T #Used in arc-length solver
    L::T #Used in dissipation solver

    #System arrays, fint, Kint, etc
    system_arrays::SystemArrays{T}

    #
    partstates::Vector{AbstractPartState}

    step::Int
    step_tries::Int
    newton_itr::Int

    #Solver specific states
    detK::T
    converged::Bool
    norm_residual::T
    solvermode::SolverMode

    #Need these for explicit solver...
    Wⁱ::T
    Wᵉ::T
    Wᵏ::T
end


function StateVariables(T::Type, ndofs::Int)
    dofvecs1 = [zeros(T, ndofs) for _ in 1:3]
    sa = SystemArrays(T, ndofs)
    return StateVariables(dofvecs1..., 0.0, 0.0, 0.0, 0.0,
                          sa, AbstractPartState[], 
                          0, 0, 0, 
                          Inf, true, Inf, MODE1, 
                          0.0, 0.0, 0.0)
end

getstatevariables_type(::AbstractSolver) = StateVariables
getsystemarrays_type(::AbstractSolver)   = SystemArrays

function transfer_state!(a::StateVariables, b::StateVariables)
    a.d .= b.d
    a.v .= b.v
    a.a .= b.a
    a.t  = b.t
    a.Δt  = b.Δt
    a.λ  = b.λ
    a.L  = b.L

    a.partstates .= b.partstates #Assume is bit-types 

    a.system_arrays = deepcopy(b.system_arrays) #This copies the stiffness matrix aswell.... hmmm

    a.step = b.step
    a.newton_itr = b.newton_itr
    a.detK = b.detK
    a.converged = b.converged #Not needed i think
    a.norm_residual = b.norm_residual
    a.solvermode = b.solvermode
end

function Base.fill!(sa::SystemArrays{T}, v::T) where T
    fill!(sa.fⁱ, v)
    fill!(sa.Kⁱ, v)
    fill!(sa.fᵉ, v)
    fill!(sa.Kᵉ, v)
    fill!(sa.fᴬ, v)
    sa.G[] = v
end


include("materials/materials.jl")
include("elements/elements.jl")

include("contact/contact.jl") #Needs reviving

include("outputs/output.jl")
include("outputs/dof_value.jl")
include("outputs/material_output.jl")
include("outputs/solverstats_output.jl")
include("outputs/energy_output.jl")

include("solvers/solver_utils.jl")
include("solvers/solver.jl")
include("solvers/dissipation_solver.jl")
include("solvers/local_dissipation_solver.jl")
include("solvers/newton_solver.jl")
include("solvers/arclength_solver.jl")
include("solvers/explicit_solver.jl")
include("solvers/implicit_solver.jl")

include("parts/parts.jl")
include("assembling.jl")

include("externalforce/external_forces.jl")
include("constraints/Constraints.jl")
#include("constraints/linearconstraints.jl)

include("utils/utils.jl")
include("utils/celliterator2.jl")
include("utils/juafemutils.jl")
include("utils/stresses_through_thickness.jl")
include("utils/adaptive.jl")

include("solvers/problem_builder.jl")

"""
    GlobalData{dim,T,DH<:AbstractDofHandler}

Contains all information about the problem being solved, e.g Forces, Boundary conditions, Parts

"""
mutable struct GlobalData{dim,T,DH<:Ferrite.AbstractDofHandler}
    dbc::ConstraintHandler{DH,T}

    grid::Ferrite.AbstractGrid{dim}
    dh::DH
    
    constraints::Constraints
    efh::ExternalForceHandler{T}
    contact::AbstractContactSearchAlgorithm
    
    parts::Vector{AbstractPart{dim}}

    output::Output{T}

    t0::T
    tend::T

    adaptive::Bool
end



end
