export DissipationOutput

struct DissipationOutput <: Five.AbstractOutput

end

function build_outputdata(output::DissipationOutput, set::Int, dh::MixedDofHandler)
    return output
end

function collect_output!(output::DissipationOutput, state::StateVariables, set::Int, globaldata)

    Δg = 0.0
    part = globaldata.parts[first(set)]

    @timeit "outputdiss" G = _copy_of_fepart(part, state, globaldata.dh)

    return G
end

function _copy_of_fepart(part::FEPart, state::StateVariables{T}, dh::Ferrite.AbstractDofHandler) where T

    element = part.element
    ue, Δue, due, ke, fe, coords, celldofs = (part.cache.ue, part.cache.Δue, part.cache.due, 
    part.cache.ke, part.cache.fe, part.cache.coords, part.cache.celldofs)

    Δt = state.Δt
    G = 0.0
    for (localid,cellid) in enumerate(part.cellset)
        
        partstate::get_partstate_type(part) = state.prev_partstates[cellid]

        materialstate = partstate.materialstates
        cellstate     = partstate.elementstate

        fill!(fe, 0.0)

        Ferrite.cellcoords!(coords, dh, cellid)
        Ferrite.celldofs!(celldofs, dh, cellid)
        
        Δue .= state.Δd[celldofs]
        ue .= state.d[celldofs]
        due .= state.v[celldofs]

        ge = Base.RefValue(zero(T))
        integrate_dissipation!(element, cellstate, part.material, materialstate, fe, ge, coords, Δue, ue, due, Δt)

        G += ge[]

    end

    return G

end
