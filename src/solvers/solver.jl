export solvethis

"""
    step!(solver::AbstractSolver, state::StateVariables, globaldata::GlobalData)

Given the current `state` of a system, the solver advances to the next state (mutating the input argument `state`).
The next state (and step size) is determined by the parameters of the `solver`.
"""
function step!(::AbstractSolver, state::StateVariables, globaldata) end

"""
    isdone(solver::AbstractSolver, state::StateVariables, globaldata::GlobalData)

Determines if the state has reached the termination criteria of the solver
"""
function Base.isdone(::AbstractSolver, state::StateVariables, globaldata) end


"""
    should_abort(solver::AbstractSolver, state::StateVariables, globaldata::GlobalData)

Determines if the solver should abort simulation (for example if it is not able to converge)
"""
should_abort(::AbstractSolver, state::StateVariables, globaldata)

function solvethis(solver::AbstractSolver, state::StateVariables, globaldata) 
    
    starttime = time()
    reset_timer!()
    
    dh, output, parts = (globaldata.dh, globaldata.output, globaldata.parts)

    #Output initial state
    @info "[SOLVER] Exporting initial state"
    export!(output, state, globaldata, force=true)
    #outputs!(output, state, globaldata, force=true)

    #Some solvers want the store the initial external force vector
    @info "[SOLVER] Init system arrays"
    @timeit "build system arrays" init_system_arrays!(solver, state, globaldata)

    prev_state = deepcopy(state)
    while !isdone(solver, state, globaldata)

        @info "[SOLVER] Step: $(state.step), Time: $(state.t)"

        success = false
        ntries = 0
        while true
            state.step += 1
            @info "[SOLVER] Calling step! $(typeof(solver)), tries: $(ntries)"
            success = step!(solver, state, globaldata, ntries)
            ntries += 1
            
            (success || should_abort(solver, state, globaldata)) && break
            
            #Reset state to previous step...
            transfer_state!(state, prev_state)
        end

        if !success
            @warn("Could not take step")
            break
        end

        @timeit "post"       post_stuff!(dh, state, globaldata)  
        @timeit "export"     export!(output, state, globaldata)
        #@timeit "vtk export" vtk_add_state!(output, state, globaldata)
        #@timeit "output"     outputs!(output, state, globaldata)

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
                Ferrite.copy!!(state.system_arrays.q, state.system_arrays.fᵉ)

                #Recalculate fᴬ for fstar
                #assemble_fstar!(globaldata.dh, state, globaldata)

                #Restart state
                prev_state = deepcopy(state)
                #continue
            end

        end
        
        #Update counter
        transfer_state!(prev_state, state)
        
        #Checkinput
        @timeit "Handle IO"  handle_input_interaction(output)
        if get_simulation_termination(output) == ABORTED_SIMULATION
            @warn("Simulation aborted by user")
            break
        end
    end

    close(LOG.stream)

    totaltime = time() - starttime
    
    print_timer(linechars = :ascii)

    @info "[SOLVER] Saving outputs"
    save_outputs(output)
    set_simulation_termination!(output, NORMAL_TERMINATION)
    
    return output
end

function reset!(a::StateVariables, b::StateVariables)
    a.d .= b.d
    a.v .= b.v
    a.a .= b.a
    a.t  = b.t
    a.λ  = b.λ
    a.L  = b.L

    a.Δd .= b.Δd
    a.Δv .= b.Δv
    a.Δa .= b.Δa
    a.Δt  = b.Δt
    a.Δλ  = b.Δλ
    a.ΔL  = b.ΔL

    a.partstates .= deepcopy(b.partstates)
    a.prev_partstates .= deepcopy(b.prev_partstates)
    
    a.step = b.step

    a.system_arrays = deepcopy(b.system_arrays)

    #Solver specific states
    a.detK = b.detK
    a.prev_detK = b.prev_detK
    a.step_tries = b.step_tries
    a.converged = b.converged
    a.norm_residual = b.norm_residual
    a.newton_itr = b.newton_itr
    a.solvermode = b.solvermode


end

function init_system_arrays!(solver::AbstractSolver, state, globaldata)
   
    fill!(state.system_arrays, 0.0)

    #Create sparsity pattern
    state.system_arrays.Kⁱ = create_sparsity_pattern(globaldata.dh)#, globaldata.dbc)
    fill!(state.system_arrays.Kⁱ.nzval, 1.0)
    for part in globaldata.parts
        K2 = assemble_sparsity_pattern!(part, globaldata)
        state.system_arrays.Kⁱ .+= K2
    end
    Ferrite._condense_sparsity_pattern!(state.system_arrays.Kⁱ, globaldata.dbc.dofcoefficients, globaldata.dbc.dofmapping, true)
    fill!(state.system_arrays.Kⁱ.nzval, 0.0)

    #TODO: Activate for time dependent solvers
    if false# !(solver isa LocalDissipationSolver) #if solver is not a ExplicitSolver
        assemble_massmatrix!(globaldata.dh, state, globaldata)
        assemble_lumped_massmatrix!(globaldata.dh, state, globaldata)
    end

    apply_external_forces!(globaldata.dh, globaldata.efh, state, globaldata)
    assemble_stiffnessmatrix_and_forcevector!(globaldata.dh, state, globaldata)
    apply_constraints!(globaldata.dh, globaldata.constraints, state,  globaldata)
    
    state.system_arrays.q .= state.system_arrays.fᵉ
    #state.prev_detK = state.detK = det(state.system_arrays.Kⁱ - state.system_arrays.Kᵉ)
end