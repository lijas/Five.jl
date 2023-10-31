module Five

using Reexport
using ForwardDiff
import LinearSolve
import Preconditioners
import IncompleteLU

import ChunkSplitters
using Parameters
using TimerOutputs
using StaticArrays
using LinearAlgebra
using SparseArrays
import Base: isdone

import Logging
import ProgressMeter
using InteractiveUtils

#import JLD2
import FileIO: save
import JLD2


@reexport using Ferrite
@reexport using Tensors
@reexport using MaterialModels

LOG = Logging.global_logger()

#function _getT() end
#function _getdim() end

const Optional{T} = Union{T,Nothing}

abstract type AbstractPart{dim} end
abstract type AbstractPartState end

@enum SolverMode MODE1 MODE2 MODE3
@enum SimulationTermination NORMAL_TERMINATION ERROR_TERMINATION ABORTED_SIMULATION _NO_TERMINATION

abstract type AbstractSolver{T} end
abstract type AbstractSystemArrays{T} end

"""
    SystemArrays(T::Type, ndofs::Int)

Contains global arrays, such as internal force vector and stiffness matrix.
"""
mutable struct SystemArrays{T} <: AbstractSystemArrays{T}
    fⁱ::Vector{T}
    Kⁱ::SparseArrays.SparseMatrixCSC{T,Int}

    fᵉ::Vector{T} #External force vector
    
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
    return SystemArrays(zeros(T,ndofs), spzeros(T,ndofs,ndofs), zeros(T,ndofs), Mᵈⁱᵃᵍ, M, zeros(T,ndofs), zeros(T,ndofs), Ref(0.0))
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
    
    d::Vector{T} # displacement
    v::Vector{T} # velocity
    a::Vector{T} # acceleration
    t::T
    Δt::T

    #Used in arc-length solvers
    λ::T 
    Δλ::T
    L::T
    ΔL::T

    #System arrays, fint, Kint, etc
    system_arrays::SystemArrays{T}

    #State variables for materials
    partstates::Vector{AbstractPartState}

    step::Int 
    step_tries::Int
    newton_itr::Int

    #Solver specific states
    detK::T #Used in crisfield solver
    converged::Bool
    norm_residual::T
    solvermode::SolverMode

    #Need these for explicit solver...
    Wⁱ::T #Internal
    Wᵉ::T #External
    Wᵏ::T #Kinetic
end


function StateVariables(T::Type, ndofs::Int)
    dofvecs1 = [zeros(T, ndofs) for _ in 1:3]
    sa = SystemArrays(T, ndofs)
    return StateVariables(dofvecs1..., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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

    a.Δλ  = b.Δλ
    a.λ  = b.λ
    a.ΔL = b.ΔL
    a.L  = b.L

    a.partstates .= deepcopy(b.partstates) #Assume is bit-types 

    a.system_arrays = deepcopy(b.system_arrays) #This copies the stiffness matrix aswell.... hmmm

    a.step = b.step
    a.newton_itr = b.newton_itr
    a.detK = b.detK
    a.converged = b.converged #Not needed i think
    a.norm_residual = b.norm_residual
    a.solvermode = b.solvermode
end

#TODO: change to fillzero!()
function zero_out_systemarrays!(sa::SystemArrays{T}) where T
    v = zero(T)
    fill!(sa.fⁱ, v)
    fill!(sa.Kⁱ, v)
    fill!(sa.fᵉ, v)
    fill!(sa.fᴬ, v)
    sa.G[] = v
end


include("materials/materials.jl")
include("materials/cohesive/matczbilinear.jl")
include("materials/cohesive/matczexponential.jl")
include("materials/phasefield_material.jl")

#Elements
include("elements/elements.jl")
include("elements/solidelement.jl")
include("elements/bar_element.jl")
include("elements/linearsolidelement.jl")
include("elements/PhaseFieldElement.jl")

#Cohesive
include("utils/cohesive_element_utils.jl")
include("elements/cohesive_element.jl")

#Contact (Not working)
include("contact/contact.jl") #Needs reviving

#Output
include("outputs/output.jl")
include("outputs/dof_value.jl")
include("outputs/material_output.jl")

#Solvers
include("solvers/solver_utils.jl")
include("solvers/solver.jl")
include("solvers/local_dissipation_solver.jl")
include("solvers/newton_solver.jl")
include("solvers/arclength_solver.jl")
# include("solvers/explicit_solver.jl")
# include("solvers/implicit_solver.jl")

#Parts
include("parts/parts.jl")
include("parts/fepart.jl")
include("parts/cohesive_part.jl")

#Forces
include("externalforce/external_forces.jl")

#Constraints
include("constraints/Constraints.jl")
include("constraints/follower_constraint.jl")

#Utils
include("assembling.jl")
include("utils/utils.jl")
include("utils/juafemutils.jl")
include("utils/adaptive.jl")
include("solvers/problem_builder.jl")


"""
    GlobalData{dim,T,DH<:AbstractDofHandler}

Contains all information about the problem being solved, e.g Forces, Boundary conditions, Parts

"""
mutable struct GlobalData{T,G<:Ferrite.AbstractGrid,DH<:Ferrite.DofHandler,CH<:Ferrite.ConstraintHandler}
    
    grid::G
    dh::DH
    dbc::CH
    
    constraints::Constraints
    efh::ExternalForceHandler
    contact::AbstractContactSearchAlgorithm
    
    parts::Vector{AbstractPart}
    output::Output

    t0::T
    tend::T
    adaptive::Bool #Not really used anymore
end

#Exports

#problem_builder.jl
export 
ProblemData,
build_problem,

#outputs
StressOutput,

#parts
Part

end
