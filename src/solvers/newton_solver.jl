export NewtonSolver

@with_kw struct NewtonSolver{T} <: AbstractSolver{T} 
    tol::T = 1e-8
    Δt0::T
    Δt_min::T = Δt0/1000
    Δt_max::T = Δt0*1000

    maxitr::Int = 10
    optitr::Int = 5
    max_residual::T = Inf
    maxitr_first_step::Int = maxitr
end

function Base.isdone(solver::NewtonSolver, state::StateVariables, globaldata)
    return state.t >= globaldata.tend
end

function should_abort(solver::NewtonSolver, state)
    return state.Δt <= solver.Δt_min #|| (state > solver.maxitr)
end

function step!(solver::NewtonSolver, state, globaldata)

    state0 = deepcopy(state)

    conv_failed = true
    ntries = 0
    Δg = 0.0
    local r
    while conv_failed
        
        set_initial_guess(solver, state, ntries)
        ntries += 1

        println("Newton solver, ntries: $(ntries), prev_state.t = $(state0.t), Δt = $(state.Δt), diss: $(state.L)")

        #Apply boundary conditions
        update!(globaldata.dbc, state.t)
        apply!(state.d, globaldata.dbc)
        state.Δd = state.d - state0.d

        state.newton_itr = 0
        while true
            state.newton_itr +=1;

            fill!(state.system_arrays, 0.0)

            @timeit "Dissipation" assemble_dissipation!(globaldata.dh, state, globaldata)
            Δg = state.system_arrays.G[]
            #Δg = 1/2 * dot(state.Δd, state0.system_arrays.q - state0.system_arrays.fᴬ)

            #Get internal force                                                                       
            @timeit "Assembling" assemble_stiffnessmatrix_and_forcevector!(globaldata.dh, state, globaldata)
            
            @timeit "ExternalForces" apply_external_forces!(globaldata.dh, globaldata.efh, state, globaldata)
            @timeit "Apply constraint" apply_constraints!(globaldata.dh, globaldata.constraints, state, globaldata)

            r = state.system_arrays.fⁱ - state.system_arrays.fᵉ
            K = state.system_arrays.Kⁱ - state.system_arrays.Kᵉ

            state.norm_residual = norm(r[Ferrite.free_dofs(globaldata.dbc)])

            #Solve 
            apply_zero!(K, r, globaldata.dbc)
            ΔΔd = K\-r

            state.Δd .+= ΔΔd
            state.d  .+= ΔΔd

            println("---->Normg: $(state.norm_residual), Δg = $(Δg)")
        
            if state.norm_residual < solver.tol
                conv_failed = false
                break
            end

            maxitr = (state.step == 1) ? (solver.maxitr_first_step) : solver.maxitr
            if state.newton_itr >= maxitr || state.norm_residual > solver.max_residual
                conv_failed = true
                break
            end

            state.partstates .= deepcopy(state0.partstates)
        end

        if conv_failed 
            if should_abort(solver, state)
                return false
            else
                #Reset the state
                copy!(state, state0)
            end
        else
            break
        end

    end
    
    state.L += Δg

    #@show norm(state0.system_arrays.fᴬ)
    #state0.system_arrays.fᴬ .= 0.0
#    assemble_fstar!(globaldata.dh, state, globaldata) 
    #state.system_arrays.fᴬ .= state0.system_arrays.fᴬ
    state.system_arrays.q .= r
    #@show norm(state.system_arrays.q)

    return true
end

function set_initial_guess(solver::NewtonSolver, state::StateVariables, ntries::Int)

    #Initilize Δt and newton_itr depending on if 
    # this is the first step
    ⁿΔt, newton_itr = (state.step == 1 || state.step == 2) ? (solver.Δt0, solver.optitr) : (state.Δt, state.newton_itr)

    if ntries == 0
        #For the first try, increase/decrease the step 
        # based on the number of newton_iteration in the previeus 
        # converged solution
        Δt = ⁿΔt * (0.5^(0.25*(newton_itr-solver.optitr)))
    else
        #If the previous newton loop failed (ie ntries != 0),
        # half the step size
        Δt = ⁿΔt * (1/2)^(ntries-1)
    end

    #Check if the timestep has been 
    # gone above or below max/min allowed timesteps
    if Δt > solver.Δt_max
        Δt = solver.Δt_max
    elseif Δt < solver.Δt_min
        Δt = solver.Δt_min
    end

    factor = Δt / ⁿΔt

    state.Δt = Δt
    state.t += Δt
    state.d += state.Δd*factor
end