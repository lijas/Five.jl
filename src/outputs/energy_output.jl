export EnergyOutput

struct EnergyOutput <: AbstractOutput

end

function build_outputdata(output::EnergyOutput, set::Nothing, dh::MixedDofHandler)
    return EnergyOutput()
end

function collect_output!(output::EnergyOutput, state::StateVariables, set::Nothing, globaldata)

    return (
        Wⁱ = state.Wⁱ,
        Wᵉ = state.Wᵉ,
        Wᵏ = state.Wᵏ,
    )
end

