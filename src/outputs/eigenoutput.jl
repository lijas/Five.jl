export EigenOutput

struct EigenOutput <: AbstractOutput

end

function build_outputdata(output::EigenOutput, set::Set{<:Any}, dh::MixedDofHandler)
    @assert isempty(set)
    return output
end

function collect_output!(output::EigenOutput, state::StateVariables, set::Set{<:Any}, globaldata)
    Kₜ = state.system_arrays.Kⁱ - state.system_arrays.Kᵉ
    apply_zero!(Kₜ, zeros(ndofs(globaldata.dh)), globaldata.dbc)
    
    vals, vecs = Arpack.eigs(Kₜ, nev=10, which = :SM)

    return (d = state.d, eval = real.(vals), evec = real.(vecs))
end
