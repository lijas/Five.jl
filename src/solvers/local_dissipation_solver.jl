export LocalDissipationSolver

const INTERNAL_ENERGY_LOCAL = MODE1
const DISSIPATION_LOCAL = MODE2

@with_kw struct LocalDissipationSolver{T} <: AbstractSolver{T}
    Δλ0::T
    ΔU_min::T
    ΔU_max::T
    λ_max::T
    λ_min::T
    ΔL0::T
    ΔL_min::T = ΔL0/10
    ΔL_max::T = ΔL0*10
    maxsteps::Int = 100
    a::T #Switch factor
    
    #Newton
    tol::T = 1.0e-9
    max_residual::T = 1e12
    optitr::Int = 5
    maxitr::Int = 10
    maxitr_first_step = maxitr
    finish_criterion::Function = finish_criterion
end

function Base.isdone(solver::LocalDissipationSolver, state::StateVariables, globaldata)
    return solver.finish_criterion(solver, state)
end

function should_abort(solver::LocalDissipationSolver, ΔP, solvermode)
    if solvermode == INTERNAL_ENERGY_LOCAL && ΔP <= solver.ΔU_min 
        return true
    elseif solvermode == DISSIPATION_LOCAL && ΔP <= solver.ΔL_min 
        return true
    end
    return false
end

function step!(solver::LocalDissipationSolver, state::StateVariables, globaldata)

    λ0     = state.λ
    q      = state.system_arrays.q
    state0 = deepcopy(state)

    #The first step we should the λ defined by the user, so take a standard Newton Step
    if state.step == 1
        success = do_first_newton_step!(solver, state, globaldata)
        return success
    end

    converged_failed = true
    ntries = 0
    


    #ΔP represents either Δτᵁ or Δλᴰ depending 
    # on the solver mode (INTERNAL_ENERGY or DISSIPATION).
    # Set ΔP to the privious time steps solution
    ΔP0 = ΔP = state.ΔL

    Δg = 0.0
    
    while converged_failed 
        ΔP = set_initial_guess!(solver, state, ΔP, ΔP0, ntries)
        
        println("Mode: $(state.solvermode), ntries: $(ntries), prev_state.λ = $(state.λ-state.Δλ), Δλ = $(state.Δλ), ΔL = $(state.ΔL)")
        
        state.newton_itr = 0
        ntries += 1
        while true
            state.newton_itr += 1
            fill!(state.system_arrays, 0.0)
            
            @timeit "Calculate dissipation" assemble_dissipation!(globaldata.dh, state, globaldata)
            state.ΔD = state.system_arrays.G[]
            state.ΔU = 0.5*dot(state.Δλ*state0.d + state0.λ*state.Δd, q)

            #Get internal force                                                                       
            @timeit "Assembling" assemble_stiffnessmatrix_and_forcevector!(globaldata.dh, state, globaldata)
            @timeit "Apply constraint" apply_constraints!(globaldata.dh, globaldata.constraints, state, globaldata)
            
            local ΔΔd, ΔΔλ, _success
            @timeit "Solve system" if state.solvermode == INTERNAL_ENERGY_LOCAL
                w = 0.5*dot(state0.d, q)
                h = 0.5*state0.λ*q

                Δg = state.ΔU
            elseif state.solvermode == DISSIPATION_LOCAL
                w = 0.0;
                h = state.system_arrays.fᴬ

                Δg = state.ΔD
            end

            #Normal stiffness matrix
            Kₜ = state.system_arrays.Kⁱ - state.system_arrays.Kᵉ
            rₜ = state.λ*q + state.system_arrays.fᵉ - state.system_arrays.fⁱ
            
            apply_zero!(Kₜ, rₜ, globaldata.dbc)

            ΔΔd, ΔΔλ, _success = _solve_dissipation_system(solver, Kₜ, rₜ, q, h, w, Δg, state.ΔL)
 
            if !_success
                converged_failed = true
                break
            end

            state.Δd += ΔΔd
            state.Δλ += ΔΔλ
            state.d  += ΔΔd
            state.λ  += ΔΔλ

            #Check convergance
            if norm(state.λ*q) <= 1e-10
                state.norm_residual = norm(rₜ[Ferrite.free_dofs(globaldata.dbc)])
            else
                state.norm_residual = norm(rₜ[Ferrite.free_dofs(globaldata.dbc)])/norm(state.λ*q)
            end
            println("------>Newton $(state.newton_itr): $(rpad("normr: $(state.norm_residual),", 20)) $(rpad("ΔU=$(state.ΔU),", 20)) $(rpad("ΔD=$(state.ΔD),", 20)) $(rpad("Δλ=$(state.Δλ),", 20)) $(rpad("ΔD/ΔU=$(state.ΔD/state.ΔU0),", 20))")

            if newton_done(solver, state.norm_residual, Δg, state.ΔL, state.solvermode)
                converged_failed = false
                break
            end

            maxitr = (state.step == 1) ? (solver.maxitr_first_step) : solver.maxitr
            if state.newton_itr >= maxitr || state.norm_residual > solver.max_residual || (abs(state.Δλ)/(solver.λ_max - solver.λ_min)) > 0.2 
                converged_failed = true
                break
            end

            #Reset partstates
            state.partstates .= deepcopy(state.prev_partstates)
        end

        if converged_failed 
            if should_abort(solver, ΔP, state.solvermode)
                return false
            else
                #Reset the state
                copy!(state, state0)
            end
        else
            break
        end
    end

    determine_solvermode!(solver, state)

    #solvermode changed in above function
    if state.solvermode == INTERNAL_ENERGY_LOCAL
        Δg = state.ΔU
    elseif state.solvermode == DISSIPATION_LOCAL
        Δg = state.ΔD
    end

    #Update the dissipation increaments for this timestep
    state.ΔL = Δg
    state.L += Δg
    state.Δt = Δg
    state.t += Δg

    return true
end

function newton_done(solver::LocalDissipationSolver, residual, Δg, ΔL, solvermode) 
    if residual < solver.tol 
        if norm(Δg-ΔL) < solver.tol
            return true
        else
            return false
        end
    end
    return false
end

function _solve_dissipation_system(solver::LocalDissipationSolver, Kₜ, rₜ, q, h, w::Float64, Δg::Float64, ΔL::Float64)

    KK = vcat(hcat(Kₜ, -q), hcat(h', w))
    ff = vcat(rₜ, -(Δg - ΔL))

    local aa
    try
        aa = KK\ff
    catch
        return 0.0 .* similar(h), 0.0, false
    end
    ΔΔd, ΔΔλ = (aa[1:end-1], aa[end])

    return ΔΔd, ΔΔλ, true
end

function set_initial_guess!(solver::LocalDissipationSolver, state::StateVariables, ΔP, ⁿΔP, ntries)

    local ΔP_max, ΔP_min

    if state.solvermode == DISSIPATION_LOCAL
        ΔP_max = solver.ΔL_max
        ΔP_min = solver.ΔL_min
    elseif state.solvermode == INTERNAL_ENERGY_LOCAL
        ΔP_max = solver.ΔU_max
        ΔP_min = solver.ΔU_min
    end

        if ntries == 0 
            #For the first try, increase/decrease the step 
            # based on the number of newton_iteration in the previeus 
            # converged solution
            ΔP *= (0.5^(0.1*(state.newton_itr-solver.optitr)))
        else
            #If the previous newton loop failed (ie ntries != 0),
            # half the step size
            ΔP /= 1.3
        end

    if ΔP < ΔP_min
        ΔP = ΔP_min
    elseif ΔP > ΔP_max
        ΔP = ΔP_max
    end

    factor = ΔP / ⁿΔP

    if state.solvermode == DISSIPATION_LOCAL
        state.ΔL = ΔP
        state.Δλ *= factor
        state.Δd *= factor 
    elseif state.solvermode == INTERNAL_ENERGY_LOCAL
        state.ΔL = ΔP
        state.Δλ *= factor
        state.Δd *= factor
    end

    state.λ += state.Δλ
    state.d += state.Δd 
    
    return ΔP
end

function determine_solvermode!(solver::LocalDissipationSolver, state::StateVariables)

    if state.solvermode == INTERNAL_ENERGY_LOCAL
        if state.ΔD > solver.a * state.ΔU

            println("Dissipation switch ΔD = $(state.ΔD) and Δu = $(state.ΔU) ")
            state.solvermode = DISSIPATION_LOCAL
            state.ΔU_negative = false;

        end
    elseif state.solvermode == DISSIPATION_LOCAL
        if state.ΔD < solver.a * state.ΔU


            println("Dissipation switch ΔD = $(state.ΔD) and Δu = $(state.ΔU) ")
            state.solvermode = INTERNAL_ENERGY_LOCAL
            state.ΔU = state.ΔU0

        end
    end

end


function do_first_newton_step!(solver, state, globaldata)

    #newton = NewtonSolver(
    #    λ0 = solver.λ0
    #)
    q      = state.system_arrays.q
    state0 = deepcopy(state)

    conv_failed = true
    
    ΔU = 0.0

    while true
        state.newton_itr +=1;

        fill!(state.system_arrays, 0.0)

        #Get internal force                                                                       
        @timeit "Assembling" assemble_stiffnessmatrix_and_forcevector!(globaldata.dh, state, globaldata)
        
        @timeit "ExternalForces" apply_external_forces!(globaldata.dh, globaldata.efh, state, globaldata)
        @timeit "Apply constraint" apply_constraints!(globaldata.dh, globaldata.constraints, state, globaldata)

        r = state.system_arrays.fⁱ - state.system_arrays.fᵉ - solver.Δλ0*state.system_arrays.q
        K = state.system_arrays.Kⁱ - state.system_arrays.Kᵉ

        state.norm_residual = norm(r[Ferrite.free_dofs(globaldata.dbc)])

        #Solve 
        apply_zero!(K, r, globaldata.dbc)
        ΔΔd = K\-r

        state.Δd .+= ΔΔd
        state.d  .+= ΔΔd
    
        println("---->Normg: $(state.norm_residual)")

        if state.norm_residual < solver.tol
            conv_failed = false
            break
        end

        if state.newton_itr >= solver.maxitr_first_step
            conv_failed = true
            break
        end

        state.partstates .= deepcopy(state0.partstates)
    end
    #Reset value for optitr
    state.newton_itr = solver.optitr

    #Record the elastic energy required for the first step
    #This elastic energy will be used to define the step length in the next step...
    state.Δλ = solver.Δλ0
    state.λ += solver.Δλ0

    ΔU = 0.5*dot(state.Δλ*state.d, q)

    if ΔU == 0
        error("Internal energy is equal to zero after the first timestep")
    end

    state.ΔU0 = ΔU
    state.ΔL = ΔU
    state.L += ΔU

    println("Internal energy: $(ΔU)")

    return !conv_failed
end