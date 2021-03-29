

"""
SolverStepStatistic

Holds some data relevant for the solvers.
"""
struct SolverStepStatistic
    λ::Float64
    L::Float64
    t::Float64
    Δλ::Float64
    ΔL::Float64
    Δt::Float64
    ntries::Int
    n_newton_itr::Int
    success::Bool
end



"""

"""
struct SolverStatOutput <: AbstractOutput
    data::Vector{SolverStepStatistic}
end

SolverStatOutput() = SolverStatOutput(SolverStepStatistic[])

function build_outputdata(output::SolverStatOutput, dh::MixedDofHandler)
    return output
end

function collect_output!(o::SolverStatOutput, state::StateVariables, globaldata)

    push!(o.data, SolverStepStatistic(state.λ, state.L, state.t, 
                                state.Δλ, state.ΔL, state.Δt, 
                                state.step_tries, state.newton_itr, state.converged))

end

should_output(o::SolverStatOutput, t) =  true #Output each step
set_last_output!(o::SolverStatOutput, t) = nothing