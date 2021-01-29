export solvethis

function solvethis(solver::AbstractSolver{T}, state::StateVariables, globaldata) where {T}
    
    starttime = time()
    reset_timer!()
    
    dh, output, parts = (globaldata.dh, globaldata.output, globaldata.parts)

    #Output initial state
    vtk_add_state!(output, state, globaldata)
    outputs!(output, state, globaldata)

    #Some solvers want the store the initial external force vector
    init_system_arrays!(solver, state, globaldata)

    ⁿstate = deepcopy(state)
    while !isdone(solver, state, globaldata)
        state.step += 1

        println("**Step $(state.step)**")
        success = step!(solver, state, globaldata)
        
        if !success
            @warn("NO SUCCESS")
            set_simulation_termination!(output, ERROR_TERMINATION)
            break
        else
            @timeit "post" post_stuff!(dh, state, globaldata)
            
            @timeit "vtk export" vtk_add_state!(output, state, globaldata)
            @timeit "output" outputs!(output, state, globaldata)

            #Currently dont have a system for adaptivity,
            # so hack the adaptivity stuff in here:
            if globaldata.adaptive
                
                @timeit "Commiting part" instructions = commit_stuff!(dh, state, globaldata)

                if length(instructions) != 0
                    #Upgrade state and prev_state
                    @timeit "Update dofhandler" update_dofhandler!(dh, state, instructions)

                    #Upgrade dirichlet conditions
                    resize!(globaldata.dbc.free_dofs, ndofs(dh)); 
                    globaldata.dbc.free_dofs .= 1:ndofs(dh)    

                    #Update q for dissipation solver
                    fill!(state.system_arrays, 0.0)
                    apply_external_forces!(dh, globaldata.efh, state, globaldata)
                    JuAFEM.copy!!(state.system_arrays.q, state.system_arrays.fᵉ)

                    #Recalculate fᴬ for fstar
                    #assemble_fstar!(globaldata.dh, state, globaldata)

                    #Restart state
                    #state = deepcopy(prev_state)
                    #continue
                end
            end

            #Update counter
            state.prev_partstates .= deepcopy(state.partstates)
            ⁿstate = deepcopy(state)
            prev_state = deepcopy(state)

            #Checkinput
            @timeit "Handle IO" handle_input_interaction(output)
            if get_simulation_termination(output) == ABORTED_SIMULATION
                @warn("Simulation aborted by user")
                break
            end
        end
    end

    close(LOG.stream)

    totaltime = time() - starttime

    print_timer(linechars = :ascii)

    save_outputs(output)
    set_simulation_termination!(output, NORMAL_TERMINATION)
    
    return output
end

function init_system_arrays!(solver::AbstractSolver, state, globaldata)
   
    fill!(state.system_arrays, 0.0)

    apply_external_forces!(globaldata.dh, globaldata.efh, state, globaldata)
    assemble_stiffnessmatrix_and_forcevector!(globaldata.dh, state, globaldata)
    apply_constraints!(globaldata.dh, globaldata.constraints, state,  globaldata)
    
    state.system_arrays.q .= state.system_arrays.fᵉ
    #state.prev_detK = state.detK = det(state.system_arrays.Kⁱ - state.system_arrays.Kᵉ)
end