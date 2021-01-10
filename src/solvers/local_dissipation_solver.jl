export LocalDissipationSolver

const INCREMENT_LOCAL = MODE1
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
end

function Base.isdone(solver::LocalDissipationSolver, state::StateVariables, globaldata)
    return state.λ > solver.λ_max || state.λ < solver.λ_min || state.step > solver.maxsteps
end

function should_abort(solver::LocalDissipationSolver, ΔP, solvermode)
    if solvermode == INCREMENT_LOCAL && ΔP <= solver.Δλ_min 
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

    #ΔP prepresents either ΔL or Δλ depending 
    # on the solver mode (INCREMENT_LOCAL or DISSIPATION_LOCAL).
    # Set ΔP to the privious time steps solution
    ΔP0 = ΔP = (state.solvermode == DISSIPATION_LOCAL ? (state.ΔL) : (state.Δλ))
    
    converged_failed = true
    ntries = 0
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

            #Get internal force                                                                       
            @timeit "Assembling" assemble_stiffnessmatrix_and_forcevector!(globaldata.dh, state, globaldata)
            @timeit "Apply constraint" apply_constraints!(globaldata.dh, globaldata.constraints, state, globaldata)

            #Normal stiffness matrix
            Kₜ = state.system_arrays.Kⁱ - state.system_arrays.Kᵉ
            rₜ = state.λ*q + state.system_arrays.fᵉ - state.system_arrays.fⁱ
            Δg = state.system_arrays.G[]
            
            apply_zero!(Kₜ, rₜ, globaldata.dbc)

            ΔΔd, ΔΔλ = _solve_dissipation_system(solver, Kₜ, rₜ, q, state.system_arrays.fᵉ, state.system_arrays.fᴬ, Δg, λ0, state.ΔL, state.solvermode)
            
            state.Δd += ΔΔd
            state.Δλ += ΔΔλ
            state.d  += ΔΔd
            state.λ  += ΔΔλ

            #Check convergance
            if norm(state.λ*q) <= 1e-10
                state.norm_residual = norm(rₜ[JuAFEM.free_dofs(globaldata.dbc)])
            else
                state.norm_residual = norm(rₜ[JuAFEM.free_dofs(globaldata.dbc)])/norm(state.λ*q)
            end
            println("------>Newton $(state.newton_itr): $(rpad("normr: $(state.norm_residual),", 30)) $(rpad("Δg=$(Δg),", 30)) $(rpad("Δλ=$(state.Δλ),", 30)) $(rpad("λ=$(state.λ),", 30))  $(rpad("maxd=$(maximum(abs.(state.d))),", 30)) $(rpad("maxd=$(maximum(abs.(state.Δd))),", 30))  $(rpad("L=$(norm(state.L)),", 30))")
            
            if newton_done(solver, state.norm_residual, Δg, state.ΔL, state.solvermode)
                converged_failed = false
                break
            end

            if state.newton_itr >= solver.maxitr || state.norm_residual > solver.max_residual
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
    state.ΔL = Δg
    state.L += Δg
    state.Δt = Δg
    state.t += Δg

    determine_solvermode!(solver, state)

    return true
end

function newton_done(solver::LocalDissipationSolver, residual, Δg, ΔL, solvermode) 
    if residual < solver.tol 
        if solvermode == DISSIPATION_LOCAL
            if norm(Δg-ΔL) < solver.tol
                return true
            else
                return false
            end
        end
        return true
    end
    return false
end

function _solve_dissipation_system(solver::LocalDissipationSolver, Kₜ, rₜ, q, fᵉ, fᴬ, Δg, λ0, ΔL, solvermode)
    local ΔΔd, ΔΔλ
    if solvermode == DISSIPATION_LOCAL
        w = 0.0
        h = fᴬ

        KK = vcat(hcat(Kₜ, -q), hcat(h', w))
        ff = vcat(rₜ, -(Δg - ΔL))

        aa = KK\ff
        ΔΔd, ΔΔλ = (aa[1:end-1], aa[end])
    elseif solvermode == INCREMENT_LOCAL
        ΔΔd  = Kₜ\rₜ
        ΔΔλ = 0.0
    else
        error("No mode")
    end

    return ΔΔd, ΔΔλ
end

function set_initial_guess!(solver::LocalDissipationSolver, state::StateVariables, ΔP, ⁿΔP, ntries)

    local ΔP_max, ΔP_min

    if state.solvermode == DISSIPATION_LOCAL
        ΔP_max = solver.ΔL_max
        ΔP_min = solver.ΔL_min
    elseif state.solvermode == INCREMENT_LOCAL
        ΔP_max = solver.Δλ_max
        ΔP_min = solver.Δλ_min
    end

    if state.step == 1
        ⁿΔP = ΔP = solver.Δλ0 * (1/2)^ntries
    else
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
    elseif state.solvermode == INCREMENT_LOCAL
        state.Δλ = ΔP
        state.Δd *= factor
    end

    state.λ += state.Δλ
    state.d += state.Δd 
    
    return ΔP
end

function determine_solvermode!(solver::LocalDissipationSolver, state::StateVariables)

    #Note: first step ΔL/abs(Δλ) == NaN, and will choose Incremenet
    #=if state.ΔL/abs(state.Δλ) < 1e-4
        return solver.solver_mode[] = INCREMENT_LOCAL
    else
        return solver.solver_mode[] = DISSIPATION_LOCAL
    end=#


    if state.solvermode == INCREMENT_LOCAL
        if state.ΔL >= solver.sw2d
            state.solvermode = DISSIPATION_LOCAL
            state.newton_itr = solver.optitr
        end
    elseif state.solvermode == DISSIPATION_LOCAL
        if state.ΔL/abs(state.Δλ) < solver.sw2i
            state.solvermode = INCREMENT_LOCAL
        end
    end


end