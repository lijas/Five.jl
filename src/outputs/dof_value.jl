export DofValueOutput

"""
    DofValueOutput(field::Symbol, comps::Vector{Int})

Outputs the displacements, internal forces, and external forces of the dofs in a given set
"""

struct DofValueOutput <: AbstractOutput
    field::Symbol
    dofs::Vector{Int} #TODO: rename to comp

    function DofValueOutput(; field::Symbol, dofs::AbstractVector{Int})
        return new(field, collect(dofs))
    end
    
    function DofValueOutput(field, dofs, fieldhandler)
        return new(field, dofs, fieldhandler)
    end

    fieldhandler::Ferrite.FieldHandler
end

#TODO: Rename all build_outpudata to init_
function build_outputdata(output::DofValueOutput, set::Set{<:Ferrite.BoundaryIndex}, dh)
    grid = dh.grid

    Ferrite._check_same_celltype(grid, cellid.(set))
    fh = getfieldhandler(dh, cellid(first(set)))
    return DofValueOutput(output.field, output.dofs, fh)
end

function collect_output!(output::DofValueOutput, state::StateVariables, set::Set{<:Ferrite.BoundaryIndex}, globaldata)

    dofs = Int[]
    for index in set
        _dofs = indexdofs(globaldata.dh, output.fieldhandler, index, output.field, output.dofs)
        append!(dofs, _dofs)
    end
    unique!(dofs)

    displacement = state.d[dofs]
    fint = state.system_arrays.fⁱ[dofs]
    fext = state.system_arrays.fᵉ[dofs]

    return (
        displacement = displacement,
        fint = fint, 
        fext = fext
    )
end
