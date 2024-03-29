export MaterialStateOutput

struct MaterialStateOutput{DT} <: AbstractOutput
    field::Symbol
    func::Function
end

function MaterialStateOutput(; field::Symbol, func::Function = mean, datatype::Type)
    return MaterialStateOutput{datatype}(field, func)
end

function build_outputdata(output::MaterialStateOutput, set::Set, dh)
    if  eltype(set) != Int 
        error("For MaterialStateOutput, only cellsets are allowed")
    end

    output
end

outputname(o::MaterialStateOutput) = string(o.field)
getdatatype(o::MaterialStateOutput{DT}) where DT = DT

function collect_output!(output::MaterialStateOutput, state::StateVariables, set::Set{Int}, globaldata)

    materialstates = Any[]
    for cellid in set
        cellmaterialstates = state.partstates[cellid].materialstates

        material_fields = Any[]
        for i in eachindex(cellmaterialstates)
            fieldvalue = getproperty(cellmaterialstates[i], output.field)
            push!(material_fields, fieldvalue)
        end

        push!(materialstates, output.func(material_fields))
    end

    return materialstates
end

struct StressOutput{DT} <: AbstractOutput
    func::Function
end

function StressOutput(; func::Function = mean)
    return StressOutput{SymmetricTensor{2,3,Float64,6}}(func)
end

function build_outputdata(output::StressOutput, set::Set, dh)
    if  eltype(set) != Int 
        error("For StressOutput, only cellsets are allowed")
    end
    
    output
end

outputname(o::StressOutput) = "σ"
getdatatype(o::StressOutput{DT}) where DT = DT

function collect_output!(output::StressOutput, state::StateVariables, set::Set{Int}, globaldata)

    cellstresses = Any[]
    for cellid in set
        partstate = state.partstates[cellid]
        matstate_qps = partstate.stresses
        _stress = output.func(matstate_qps)

        push!(cellstresses, _stress)
    end

    return cellstresses
end