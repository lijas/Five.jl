export MaterialStateOutput

struct MaterialStateOutput <: AbstractOutput
    field::Symbol
    func::Function
end

function MaterialStateOutput(; field::Symbol, func::Function = maximum)
    return MaterialStateOutput(field, func)
end

function collect_output!(output::MaterialStateOutput, state::StateVariables, set::Set{<:Index}, globaldata)

    for cellid in set
        partstate = state.partcellstates[cellid]
        matstate_qps = getmaterialstate(partstate, output.field)
        matstate_value = output.func(matstate_qps)
    end

    return (
        output.field = matstate_value
    )
end