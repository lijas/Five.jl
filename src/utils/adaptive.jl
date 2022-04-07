export FieldDimUpgradeInstruction

abstract type AbstractUpdateInstruction end

struct FieldDimUpgradeInstruction
    cellid::Int
    nodeids::Vector{Int}
    local_dof_idx::Vector{Int}
    current_fielddim::Vector{Int}
    upgrade_fielddim::Vector{Int}
    extended_ue::Vector{Vector{Float64}}
    extended_Δue::Vector{Vector{Float64}}
end


function update_dofhandler!(dh::MixedDofHandler, state::StateVariables{T}, instructions::Vector{FieldDimUpgradeInstruction}) where T

    vertexdict = Dict{Int, Array{Int}}()

    for instr in instructions
        _update!(dh, state, instr, vertexdict)
    end

    resize!(state.system_arrays.fᵉ, ndofs(dh))
    resize!(state.system_arrays.fⁱ, ndofs(dh))
    resize!(state.system_arrays.fᴬ, ndofs(dh))
    state.system_arrays.Kᵉ = spzeros(T, ndofs(dh), ndofs(dh))
    state.system_arrays.Kⁱ = create_sparsity_pattern(dh)

    #Resize the velocites anc accelerations since they are not used
    resize!(state.v, ndofs(dh))
    resize!(state.a, ndofs(dh))

end

function _update!(dh::MixedDofHandler, state::StateVariables, instr::FieldDimUpgradeInstruction, vertexdict)

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
        nextdof, dofs = Ferrite.get_or_create_dofs!(nextdof, Δ, dict=vertexdict, key=instr.nodeids[i])
        for j in 0:Δ-1
            insert!(dh.cell_dofs.values, celloffset-1 + instr.local_dof_idx[i] + j, dofs[j+1])
        end

        #Resize and add displacement values for the new dofs
        resize!(state.d, nextdof-1)
        state.d[nodedofs] = instr.extended_ue[i][1:current_fielddim]
        state.d[dofs]     = instr.extended_ue[i][(1:length(dofs)) .+ current_fielddim]

        resize!(state.v, nextdof-1)
        state.v[nodedofs] = instr.extended_Δue[i][1:current_fielddim]
        state.v[dofs]     = instr.extended_Δue[i][(1:length(dofs)) .+ current_fielddim]

        dh.cell_dofs.length[instr.cellid] += Δ
    end

    for idx in (instr.cellid+1):length(dh.cell_dofs.offset)
        dh.cell_dofs.offset[idx] += Δ_total
    end
    dh.ndofs[] = nextdof-1
    

end