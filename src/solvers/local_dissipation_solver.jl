export LocalDissipationSolver

const ELASTIC_LOCAL = MODE1
const DISSIPATION_LOCAL = MODE2

@with_kw struct LocalDissipationSolver{T} <: AbstractSolver{T}
    Δλ0::T
    Δλ_min = Δλ0/1000
    Δλ_max = Δλ0*1000
    λ_max::T
    λ_min::T
    ΔL0::T
    ΔL_min::T = ΔL0/10
    ΔL_max::T = ΔL0*10
    maxsteps::Int = 100
    sw2d::T
    sw2i::T
    
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
    if solvermode == ELASTIC_LOCAL && ΔP <= solver.Δλ_min 
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

    #
    if state.step == 1
        success = do_first_newton_step!(solver, state, globaldata)
        return success
    end

    converged_failed = true
    ntries = 0
    
    Δg = state.Δt # Time variable is borrowed to store the dissipation for this solver.
    determine_solvermode!(solver, state, Δg)

    #ΔP represents either ΔL or Δλ depending 
    # on the solver mode (INCREMENT or DISSIPATION).
    # Set ΔP to the privious time steps solution
    ΔP0 = ΔP = state.ΔL


    while converged_failed 
        ΔP = set_initial_guess!(solver, state, ΔP, ΔP0, ntries)
        
        println("Mode: $(state.solvermode), ntries: $(ntries), prev_state.λ = $(state.λ-state.Δλ), Δλ = $(state.Δλ), ΔL = $(state.ΔL)")
        
        state.newton_itr = 0
        ntries += 1
        while true
            state.newton_itr += 1
            fill!(state.system_arrays, 0.0)
            
            @timeit "Calculate dissipation" assemble_dissipation!(globaldata.dh, state, globaldata)
            Δg = state.system_arrays.G[]

            if state.solvermode == ELASTIC_LOCAL
                fill!(state.system_arrays, 0.0)
                assemble_elastic!(globaldata.dh, state, globaldata)
            end   
            
            #Get internal force                                                                       
            @timeit "Assembling" assemble_stiffnessmatrix_and_forcevector!(globaldata.dh, state, globaldata)
            @timeit "Apply constraint" apply_constraints!(globaldata.dh, globaldata.constraints, state, globaldata)
            
            #Normal stiffness matrix
            Kₜ = state.system_arrays.Kⁱ - state.system_arrays.Kᵉ
            rₜ = state.λ*q + state.system_arrays.fᵉ - state.system_arrays.fⁱ
            
            apply_zero!(Kₜ, rₜ, globaldata.dbc)

            @timeit "Solve system" ΔΔd, ΔΔλ, _success = _solve_dissipation_system(solver, Kₜ, rₜ, q, state.system_arrays.fᴬ, state.system_arrays.G[], state.ΔL)

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
            #println("------>Newton $(state.newton_itr): $(rpad("normr: $(state.norm_residual),", 20)) $(rpad("Δg=$(state.system_arrays.G[]),", 20)) $(rpad("Δλ=$(state.Δλ),", 20)) $(rpad("λ=$(state.λ),", 20))  $(rpad("maxd=$(maximum(abs.(state.d))),", 30)) $(rpad("maxd=$(maximum(abs.(state.Δd))),", 30))  $(rpad("L=$(norm(state.L)),", 30))")
            @printf("--->Newton %i, normr: %.5e, dissi: %.5e, Δλ: %.5e, Δg:%.5e \n", state.newton_itr, state.norm_residual, Δg, state.Δλ, state.system_arrays.G[])

            if newton_done(solver, state.norm_residual, state.system_arrays.G[], state.ΔL, state.solvermode)
                converged_failed = false
                break
            end

            maxitr = (state.step == 1) ? (solver.maxitr_first_step) : solver.maxitr
            if state.newton_itr >= maxitr || state.norm_residual > solver.max_residual || (abs(state.Δλ)/(solver.λ_max - solver.λ_min)) > 2.0 
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

    #Update the dissipation increaments for this timestep
    state.ΔL = state.system_arrays.G[]
    state.L += state.system_arrays.G[]
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

function _solve_dissipation_system(solver::LocalDissipationSolver, Kₜ, rₜ, q, fᴬ, Δg, ΔL)
    local ΔΔd, ΔΔλ

    w = 0.0
    h = fᴬ

    KK = vcat(hcat(Kₜ, -q), hcat(h', w))
    ff = vcat(rₜ, -(Δg - ΔL))
    aa = KK\ff
    local aa
    try
        aa = KK\ff
    catch
        return 0.0 .* copy(h), 0.0, false
    end
    
    ΔΔd, ΔΔλ = (aa[1:end-1], aa[end])

    return ΔΔd, ΔΔλ, true
end

function set_initial_guess!(solver::LocalDissipationSolver, state::StateVariables, ΔP, ⁿΔP, ntries)

    local ΔP_max, ΔP_min

    if state.solvermode == DISSIPATION_LOCAL
        ΔP_max = solver.ΔL_max
        ΔP_min = solver.ΔL_min
    elseif state.solvermode == ELASTIC_LOCAL
        ΔP_max = solver.Δλ_max
        ΔP_min = solver.Δλ_min
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
    
    state.ΔL = ΔP
    state.Δλ *= factor
    state.Δd *= factor
    
    state.λ += state.Δλ
    state.d += state.Δd 
    
    return ΔP
end

function determine_solvermode!(solver::LocalDissipationSolver, state::StateVariables, Δg)


    if state.solvermode == ELASTIC_LOCAL
        if Δg >= solver.sw2d
            state.solvermode = DISSIPATION
            #state.newton_itr = solver.optitr
            state.ΔL = Δg
        end
    elseif state.solvermode == DISSIPATION_LOCAL
        #error("No switching back")
        #if Δg/abs(state.Δλ) < solver.sw2i
        #    state.solvermode = ELASTIC_LOCAL
        #end
    end


end



function do_first_newton_step!(solver, state, globaldata)

    #newton = NewtonSolver(
    #    λ0 = solver.λ0
    #)

    state0 = deepcopy(state)

    conv_failed = true
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

    assemble_elastic!(globaldata.dh, state, globaldata)
    state.system_arrays.G[] = 1e-6
    
    state.newton_itr = solver.optitr

    state.ΔL = state.system_arrays.G[]
    state.L += state.system_arrays.G[]
    state.Δλ = solver.Δλ0
    state.λ = solver.Δλ0

    @show state.ΔL

    return !conv_failed
end