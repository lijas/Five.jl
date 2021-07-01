export ImplicitSolver

@with_kw struct ImplicitSolver{T} <: AbstractSolver{T} 
    Δt0::T
    β::T
    γ::T
    tol::T

    optitr = 5
    max_residual::T = Inf
    maxitr::Int = 10
end

function Base.isdone(solver::ImplicitSolver, state::StateVariables, globaldata)
    return state.t >= globaldata.tend
end

function should_abort(solver::ImplicitSolver, state)
    return false
end


function step!(solver::ImplicitSolver{T}, state, globaldata) where {T}
    
    state0 = deepcopy(state)

    
    M = state.system_arrays.M
    converged = false
    ntries = 0

    while !converged
        set_initial_guess(solver, state, ntries)
        ntries += 1
        Δt = state.Δt

        println("Implicit solver, ntries: $(ntries), prev_state.t = $(state0.t), Δt = $(state.Δt)")

        state.t += Δt

        d̃ = state.d + Δt*state.v + (Δt^2)/2*(1 - 2*solver.β)*state.a;
        ṽ = state.v + (1-solver.γ)*Δt*state.a;
    
        update!(globaldata.dbc, state.t)
        
        state.d .= d̃
        apply!(state.d, globaldata.dbc)

        state.newton_itr = 0 
        while true
            state.newton_itr +=1;
            fill!(state.system_arrays, 0.0)

            #Get internal force                                                                       
            @timeit "Assembling" assemble_stiffnessmatrix_and_forcevector!(globaldata.dh, state, globaldata)
            
            @timeit "ExternalForces" apply_external_forces!(globaldata.dh, globaldata.efh, state, globaldata)
            @timeit "Apply constraint" apply_constraints!(globaldata.dh, globaldata.constraints, state, globaldata)

            f = state.system_arrays.fᵉ - state.system_arrays.fⁱ

            #comute a_np1 and v_np1
            state.a = 1/(solver.β*Δt^2)*(state.d - d̃);
            state.v = ṽ + solver.γ * Δt * state.a;

            #Compute residule
            r = M*state.a - f;

            state.norm_residual = norm(r[free_dofs(globaldata.dbc)])

            #Compute jacobian
            A = 1/(solver.β*Δt^2)*M + state.system_arrays.Kⁱ;

            #Solve 
            apply_zero!(A, r, globaldata.dbc)
            #try
                @timeit "Solve" ΔΔd = A\-r
            #catch e
            #    println(e)
            #    conveged = false
            #end
            #apply_zero!(ΔΔd, globaldata.dbc)

            #updata
            state.Δd += ΔΔd
            state.d  += ΔΔd

            println("---->Normg: $(state.norm_residual)")

            scaledtol = solver.tol * max( norm(state.system_arrays.fⁱ), norm(state.system_arrays.fᵉ), norm(M*state.a) )
            #Check convergence
            if state.norm_residual < scaledtol
                converged = true
                break
            end

            if state.newton_itr >= solver.maxitr || state.norm_residual > solver.max_residual
                converged = false
                break
            end

            state.partstates .= deepcopy(state0.partstates)
        end

        if !converged 
            if should_abort(solver, state)
                return false
            else
                #Reset the state
                copy!(state, state0)
                continue
            end
        else
            break
        end

    end
      
    return converged
end


function set_initial_guess(solver::ImplicitSolver, state::StateVariables, ntries::Int)


    ⁿΔt, newton_itr = (state.step == 1) ? (solver.Δt0, solver.optitr) : (state.Δt, state.newton_itr)

    if ntries == 0
        #For the first try, increase/decrease the step 
        # based on the number of newton_iteration in the previeus 
        # converged solution
        Δt = ⁿΔt * (0.5^(0.1*(newton_itr-solver.optitr)))
    else
        #If the previous newton loop failed (ie ntries != 0),
        # half the step size
        Δt = ⁿΔt * (1/2)^(ntries)
    end

    #if Δt > solver.Δt_max
    #    Δt = solver.Δt_max
    #elseif Δt < solver.Δt_min
    #    Δt = solver.Δt_min
    #end

    #factor = Δt / ⁿΔt

    state.Δt = Δt
    #state.d += state.Δd*factor

end