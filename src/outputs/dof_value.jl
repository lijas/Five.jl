export DofValueOutput

struct DofValueOutput <: AbstractOutput
    field::Symbol
    dofs::Vector{Int}

    function DofValueOutput(; field::Symbol, dofs::AbstractVector{Int})
        return new(field, collect(dofs))
    end
    
    function DofValueOutput(field, dofs, fieldhandler)
        return new(field, dofs, fieldhandler)
    end

    fieldhandler::JuAFEM.FieldHandler
end

function build_outputdata(output::DofValueOutput, set::Set{<:Index}, dh::MixedDofHandler)
    JuAFEM._check_same_celltype(dh.grid, cellid.(set))
    fh = getfieldhandler(dh, cellid(first(set)))
    return DofValueOutput(output.field, output.dofs, fh)
end

function collect_output!(output::DofValueOutput, state::StateVariables, set::Set{<:Index}, globaldata)

    dofs = Int[]
    for index in set
        _dofs = indexdofs(globaldata.dh, output.fieldhandler, index, output.field, output.dofs)
        append!(dofs, _dofs)
    end
    unique!(dofs)

    displacement = maximum(abs.(state.d[dofs]))
    fint = sum(state.system_arrays.fⁱ[dofs])
    fext = sum(state.system_arrays.fᵉ[dofs])

    return (
        displacement = displacement,
        fint = fint, 
        fext = fext
    )
end

function collect_output!(output::DofValueOutput, state::StateVariables, set::Index, globaldata) 

    dof = first(dofs_on_vertex(dh, output.fieldhandler, output.vertexindex, output.field, [output.component]))

    displacement = maximum(abs.(state.d[dofs]))
    fint = sum(state.system_arrays.fⁱ[dofs])
    fext = sum(state.system_arrays.fᵉ[dofs])

    return (
        displacement = displacement,
        fint = fint, 
        fext = fext
    )
end
