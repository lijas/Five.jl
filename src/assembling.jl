

function assemble_massmatrix!(dh, state, globaldata)
    for (partid, part) in enumerate(globaldata.parts)
        assemble_massmatrix!(dh, part, state)
    end
end

function assemble_lumped_massmatrix!(dh, state, globaldata)

    assemble_massmatrix!(dh, system_arrays, globaldata)

    fill!(system_arrays.Mᵈⁱᵃᵍ, 0.0)
    for I in 1:ndofs(dh)
        for J in 1:ndofs(dh)
            system_arrays.Mᵈⁱᵃᵍ[I,I] += system_arrays.M[I,J]
        end
    end
end

function update_massmatrix!(dh, state, globaldata)

    for (partid, part) in enumerate(parts)
        
        if has_constant_massmatrix(part.element)
            continue
        end

        material = materials[part.materialid]
        elementinfo = elementinfos[part.elementinfo_id]

        ndofs = Ferrite.ndofs(part.element)
        me = zeros(ndofs,ndofs)

        for cell in CellIterator(dh, part.element, part.cellset)
            
            fill!(me, 0)
            elementstate = cellstates[cell.current_cellid[]]

            global_dofs = cell.celldofs

            integrate_massmatrix!(part.element, elementstate, elementinfo, material, cell, me, u[global_dofs], du[global_dofs])
            
            #This does not work generally, but it works for rigid bodies since they do not share dofs with other elements
            M[global_dofs,global_dofs] = me

        end
    end
end

function assemble_bodyforce!(dh, parts::Vector{<:AbstractPart}, state)
    for (partid, part) in enumerate(parts)
        assemble_bodyforce!(dh, part, materials, cellstates, elementinfos, f, forcevec)
    end
end

                                                  
function assemble_stiffnessmatrix_and_forcevector!(dh, state::StateVariables, globaldata)
    for (partid, part) in enumerate(globaldata.parts)
        assemble_stiffnessmatrix_and_forcevector!(dh, part, state)
    end
end

function assemble_forcevector!(dh, state::StateVariables, globaldata)
    
    for (partid, part) in enumerate(globaldata.parts)

        assemble_forcevector!(dh, part, state)
    end
end

function assemble_fstar!(dh, state, globaldata)
    
    for (partid, part) in enumerate(globaldata.parts)

        assemble_fstar!(dh, part, state)
    end
end

function assemble_dissipation!(dh::MixedDofHandler, state::StateVariables, globaldata)
    for (partid, part) in enumerate(globaldata.parts)

        assemble_dissipation!(dh, 
                                part, 
                                state)
    end
end


function post_stuff!(dh, state, globaldata)
    for (partid, part) in enumerate(globaldata.parts)
        post_part!(dh, part, state)
    end
end

#For adaptivity, future
function commit_stuff!(dh, state, globaldata)
    instructions = FieldDimUpgradeInstruction[]
    for (partid, part) in enumerate(globaldata.parts)
        instr = commit_part!(dh, part, state)
        if instr !== nothing
            append!(instructions, instr)
        end
    end
    return instructions
end
