export DofValueOutput

"""
    DofValueOutput(field::Symbol, comps::Vector{Int})

Outputs the displacements, internal forces, and external forces of the dofs in a given set
"""

struct DofValueOutput <: AbstractOutput
    field::Symbol
    dofs::Vector{Int}

    function DofValueOutput(; field::Symbol, dofs::AbstractVector{Int})
        return new(field, collect(dofs))
    end
    
    function DofValueOutput(field, dofs, subdofhandler)
        return new(field, dofs, subdofhandler)
    end
end

function build_outputdata(output::DofValueOutput, set::Set{<:Ferrite.BoundaryIndex}, dh)
    grid = dh.grid

    Ferrite._check_same_celltype(grid, cellid.(set))
    sdh = getsubdofhandler(dh, cellid(first(set)))

    for index in set
        _dofs = indexdofs(dh, sdh, index, output.field, output.dofs)
        append!(output.dofs, _dofs)
    end
    unique!(output.dofs)

    #return DofValueOutput(output.field, output.dofs)
    return output
end

function collect_output!(output::DofValueOutput, state::StateVariables, set::Set{<:Ferrite.BoundaryIndex}, globaldata)

    dofs = output.dofs

    displacement = state.d[dofs]
    fint = state.system_arrays.fⁱ[dofs]
    fext = state.system_arrays.fᵉ[dofs]

    return (
        displacement = displacement,
        fint = fint, 
        fext = fext
    )
end
