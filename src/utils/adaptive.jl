export FieldDimUpgradeInstruction

abstract type AbstractUpdateInstruction end

struct FieldDimUpgradeInstruction
    cellid::Int
    nodeids::Vector{Int}
    local_dof_idx::Vector{Int}
    current_fielddim::Vector{Int}
    upgrade_fielddim::Vector{Int}
    extended_ue::Vector{Vector{Float64}}
    extended_ue_prev::Vector{Vector{Float64}}
    extended_Δue::Vector{Vector{Float64}}
    extended_Δue_prev::Vector{Vector{Float64}}
end


function update_dofhandler!(dh::MixedDofHandler, state::StateVariables{T}, prev_state::StateVariables, systemarrays::SystemArrays, instructions::Vector{FieldDimUpgradeInstruction}) where T

    vertexdict = Dict{Int, Array{Int}}()

    for instr in instructions
        _update!(dh, state, prev_state, systemarrays, instr, vertexdict)
    end

    resize!(systemarrays.fᵉ, ndofs(dh))
    resize!(systemarrays.fⁱ, ndofs(dh))
    resize!(systemarrays.fᴬ, ndofs(dh))
    systemarrays.Kᵉ = spzeros(T, ndofs(dh), ndofs(dh))
    systemarrays.Kⁱ = create_sparsity_pattern(dh)

    #Resize the velocites anc accelerations since they are not used
    resize!(state.v, ndofs(dh))
    resize!(state.a, ndofs(dh))

    error("Update other")
end

function _update!(dh::MixedDofHandler, state::StateVariables, prev_state::StateVariables, systemarrays::SystemArrays, instr::FieldDimUpgradeInstruction, vertexdict)

    _celldofs = celldofs(dh, instr.cellid)

    #Modify celldofs vector
    nextdof = ndofs(dh) + 1
    Δ_total = 0
    offset = 0

    # Must sort local_dof_idx since we must insert dofs in reverse
    perm = sortperm(instr.local_dof_idx, rev=true)

    for i in perm
        #println("  Node $(instr.nodeids[i]), dofoffset: $(instr.local_dof_idx[i])")
        current_fielddim = instr.current_fielddim[i]

        #Get the celldofs for the current node
        nodedofs = _celldofs[(1:current_fielddim) .+ (instr.local_dof_idx[i] -1 - current_fielddim)]
        offset += current_fielddim

        #println("  Updating $current_fielddim -> $(instr.upgrade_fielddim[i])")

        #Get the number of new dofs, Δ, for this node
        Δ = instr.upgrade_fielddim[i] - current_fielddim
        Δ_total += Δ
        celloffset = dh.cell_dofs.offset[instr.cellid]

        #Get the new dofs, and insert them
        nextdof, dofs = JuAFEM.get_or_create_dofs!(nextdof, Δ, dict=vertexdict, key=instr.nodeids[i])
        for j in 0:Δ-1
            insert!(dh.cell_dofs.values, celloffset-1 + instr.local_dof_idx[i] + j, dofs[j+1])
        end

        #Resize and add displacement values for the new dofs
        resize!(state.d, nextdof-1)
        state.d[nodedofs] = instr.extended_ue[i][1:current_fielddim]
        state.d[dofs]     = instr.extended_ue[i][(1:length(dofs)) .+ current_fielddim]

        resize!(prev_state.d, nextdof-1)
        prev_state.d[nodedofs] = instr.extended_ue_prev[i][1:current_fielddim]
        prev_state.d[dofs]     = instr.extended_ue_prev[i][(1:length(dofs)) .+ current_fielddim]

        resize!(state.Δd, nextdof-1)
        state.Δd[nodedofs] = instr.extended_Δue[i][1:current_fielddim]
        state.Δd[dofs]     = instr.extended_Δue[i][(1:length(dofs)) .+ current_fielddim]

        resize!(prev_state.Δd, nextdof-1)
        prev_state.Δd[nodedofs] = instr.extended_Δue_prev[i][1:current_fielddim]
        prev_state.Δd[dofs]     = instr.extended_Δue_prev[i][(1:length(dofs)) .+ current_fielddim]

        dh.cell_dofs.length[instr.cellid] += Δ
    end

    for idx in (instr.cellid+1):length(dh.cell_dofs.offset)
        dh.cell_dofs.offset[idx] += Δ_total
    end
    dh.ndofs[] = nextdof-1
    

end



function initial_upgrade_of_dofhandler(dh::MixedDofHandler, igashell::IGAShell)

    instructions = FieldDimUpgradeInstruction[]

    for (ic, cellid) in enumerate(igashell.cellset)

        cellnodes = igashell.cell_connectivity[:, ic]
        cellnode_states = @view adapdata(igashell).control_point_states[cellnodes]

        initial_cellnode_states = fill(LUMPED, length(cellnode_states))

        if cellnode_states != initial_cellnode_states
            ndofs = ndofs_per_cell(dh, cellid)

            instr = construct_upgrade_instruction(igashell, cellid, initial_cellnode_states, cellnode_states, zeros(Float64, ndofs), zeros(Float64, ndofs), zeros(Float64, ndofs), zeros(Float64, ndofs))
            push!(instructions, instr)
        end

    end
    return instructions
        
end