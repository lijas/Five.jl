module Five

using Reexport
using ForwardDiff
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
import IGA #https://github.com/lijas/IGA.jl.git


@reexport using JuAFEM
@reexport using Tensors

LOG = Logging.global_logger()

abstract type AbstractPart{dim} end
abstract type AbstractPartState end

@enum SolverMode MODE1 MODE2
@enum SimulationTermination NORMAL_TERMINATION ERROR_TERMINATION ABORTED_SIMULATION _NO_TERMINATION

abstract type AbstractSolver{T} end

export SystemArrays, StateVariables, GlobalData

mutable struct SystemArrays{T}
    fⁱ::Vector{T}
    Kⁱ::SparseArrays.SparseMatrixCSC{T,Int}

    fᵉ::Vector{T}
    Kᵉ::SparseArrays.SparseMatrixCSC{T,Int}
    
    Mᵈⁱᵃᵍ::SparseArrays.SparseMatrixCSC{T,Int}
    M    ::SparseArrays.SparseMatrixCSC{T,Int}
    
    #For arc-solvers
    q::Vector{T}

    #For dissipation solver:
    fˢ::Vector{T}
    fᴬ::Vector{T}

    G::Base.RefValue{T}
end

function SystemArrays(T::Type, ndofs::Int)
    Mᵈⁱᵃᵍ = spzeros(T,ndofs,ndofs)
    M = spzeros(T,ndofs,ndofs)
    return SystemArrays(zeros(T,ndofs), spzeros(T,ndofs,ndofs), zeros(T,ndofs), spzeros(T,ndofs,ndofs), Mᵈⁱᵃᵍ, M, zeros(T,ndofs), zeros(T,ndofs), zeros(T,ndofs), Ref(0.0))
end

mutable struct StateVariables{T}
    
    d::Vector{T}
    v::Vector{T}
    a::Vector{T}
    t::T
    λ::T #Used in arc-length solver
    L::T #Used in dissipation solver

    Δd::Vector{T}
    Δv::Vector{T}
    Δa::Vector{T} 
    Δt::T
    Δλ::T
    ΔL::T

    #System arrays, fint, Kint, etc
    system_arrays::SystemArrays{T}
    #prev_system_arrays::SystemArrays

    #A bit difficult to store partstates and Δpartstates,
    # so store current and previous time of partstates instead. 
    partstates::Vector{AbstractPartState}
    prev_partstates::Vector{AbstractPartState}

    step::Int

    #Solver specific states
    detK::T
    prev_detK::T
    step_tries::Int
    converged::Bool
    norm_residual::T
    newton_itr::Int
    solvermode::SolverMode
end

function StateVariables(T::Type, ndofs::Int)
    dofvecs1 = [zeros(T, ndofs) for _ in 1:3]
    dofvecs2 = [zeros(T, ndofs) for _ in 1:3]
    sa = SystemArrays(T, ndofs)
    return StateVariables(dofvecs1..., 0.0, 0.0, 0.0, dofvecs2..., 0.0, 0.0, 0.0, sa, AbstractPartState[], AbstractPartState[], 0, NaN, NaN, 0, true, Inf, 0, MODE1)
end

function Base.copy!(a::StateVariables, b::StateVariables)
    a.d .= b.d
    a.v .= b.v
    a.a .= b.a
    a.t  = b.t
    a.λ  = b.λ
    a.L  = b.L

    a.Δd .= b.Δd
    a.Δv .= b.Δv
    a.Δa .= b.Δa
    a.Δt  = b.Δt
    a.Δλ  = b.Δλ
    a.ΔL  = b.ΔL

    a.partstates .= deepcopy(b.partstates)
    a.prev_partstates .= deepcopy(b.prev_partstates)
    
    a.step = b.step

    a.system_arrays = deepcopy(b.system_arrays)

    #Solver specific states
    a.detK = b.detK
    a.prev_detK = b.prev_detK
    a.step_tries = b.step_tries
    a.converged = b.converged
    a.norm_residual = b.norm_residual
    a.newton_itr = b.newton_itr
    a.solvermode = b.solvermode
end

function Base.fill!(sa::SystemArrays{T}, v::T) where T
    fill!(sa.fⁱ, v)
    fill!(sa.Kⁱ, v)
    fill!(sa.fᵉ, v)
    fill!(sa.Kᵉ, v)
end


include("materials/materials.jl")
include("elements/elements.jl")

include("contact/contact.jl") #Needs reviving

include("outputs/output.jl")
include("outputs/dof_value.jl")
include("outputs/material_output.jl")

include("solvers/solver.jl")
include("solvers/dissipation_solver.jl")
include("solvers/local_dissipation_solver.jl")
include("solvers/newton_solver.jl")
include("solvers/arclength_solver.jl")

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

mutable struct GlobalData{dim,T,DH<:JuAFEM.AbstractDofHandler}
    dbc::ConstraintHandler{DH,T}

    grid::Grid{dim}
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
