export DissipationSolver

const INCREMENT = MODE1
const DISSIPATION = MODE2

@with_kw struct DissipationSolver{T} <: AbstractSolver{T}
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

function Base.isdone(solver::DissipationSolver, state::StateVariables, globaldata)
    return state.λ > solver.λ_max || state.λ < solver.λ_min || state.step > solver.maxsteps
end

function should_abort(solver::DissipationSolver, ΔP, solvermode)
    if solvermode == INCREMENT && ΔP <= solver.Δλ_min 
        return true
    elseif solvermode == DISSIPATION && ΔP <= solver.ΔL_min 
        return true
    end
    return false
end

function step!(solver::DissipationSolver, state::StateVariables, globaldata)

    λ0     = state.λ
    q      = state.q
    state0 = deepcopy(state)

    #ΔP prepresents either ΔL or Δλ depending 
    # on the solver mode (INCREMENT or DISSIPATION).
    # Set ΔP to the privious time steps solution
    ΔP0 = ΔP = (state.solvermode == DISSIPATION ? (state.ΔL) : (state.Δλ))
    
    converged_failed = true
    ntries = 0
    Δg = 0.0
    #println("------>Before $(state.newton_itr): $(rpad("normr: $(state.norm_residual),", 30)) $(rpad("Δg=$(Δg),", 30)) $(rpad("Δλ=$(state.Δλ),", 30)) $(rpad("λ=$(state.λ),", 30))  $(rpad("maxd=$(maximum(abs.(state.d))),", 30)) $(rpad("maxd=$(maximum(abs.(state.Δd))),", 30))  $(rpad("fs=$(norm(state.fˢ)),", 30))")
    while converged_failed 
        ΔP = set_initial_guess!(solver, state, ΔP, ΔP0, ntries)
        
        println("Mode: $(state.solvermode), ntries: $(ntries),  Δλ = $(state.Δλ), ΔL = $(state.ΔL)")
        
        state.newton_itr = 0
        norm_residual_prev = Inf
        ntries += 1
        while true
            state.newton_itr += 1
            fill!(system_arrays, 0.0)
            #Get internal force                                                                       
            assemble_stiffnessmatrix_and_forcevector!(globaldata.dh, state, globaldata)
            apply_constraints!(globaldata.dh, globaldata.constraints, state,  globaldata)
            
            #Normal stiffness matrix
            Kₜ = state.system_arrays.Kⁱ - state.system_arrays.Kᵉ
            rₜ = state.λ*q + state.system_arrays.fᵉ - state.system_arrays.fⁱ
            
            Δg = 1/2 * dot(state.Δd, λ0*q - state.fˢ)
            
            #println("-----|>Newton $(state.newton_itr): $(rpad("normr: $(state.norm_residual),", 30)) $(rpad("Δg=$(Δg),", 30)) $(rpad("Δλ=$(state.Δλ),", 30)) $(rpad("λ=$(state.λ),", 30))  $(rpad("maxd=$(maximum(abs.(state.d))),", 30)) $(rpad("maxd=$(maximum(abs.(state.Δd))),", 30))  $(rpad("fs=$(norm(state.fˢ)),", 30))")
            apply_zero!(Kₜ, rₜ, globaldata.dbc)

            ΔΔd, ΔΔλ = _solve_dissipation_system(Kₜ, rₜ, q, state.system_arrays.fᵉ, state.fˢ, Δg, λ0, state.ΔL, state.solvermode)
            
            state.Δd += ΔΔd
            state.Δλ += ΔΔλ
            state.d  += ΔΔd
            state.λ  += ΔΔλ

            #Arc-leanth equation
            Δg = 1/2 * dot(state.Δd, λ0*q - state.fˢ)# - state.ΔL
            #Check convergance
            state.norm_residual = norm(rₜ[JuAFEM.free_dofs(globaldata.dbc)])#/norm(state.λ*state.q)
            println("------>Newton $(state.newton_itr): $(rpad("normr: $(state.norm_residual),", 30)) $(rpad("Δg=$(Δg),", 30)) $(rpad("Δλ=$(state.Δλ),", 30)) $(rpad("λ=$(state.λ),", 30))  $(rpad("maxd=$(maximum(abs.(state.d))),", 30)) $(rpad("maxd=$(maximum(abs.(state.Δd))),", 30))  $(rpad("L=$(norm(state.L)),", 30))")
            #println("Maximum: $(maximum(abs.(state.d)))")
            if newton_done(state.norm_residual, Δg, state.ΔL, solver.tol, state.solvermode)
                converged_failed = false
                break
            end

            if (state.newton_itr >= solver.maxitr || state.norm_residual > solver.max_residual) && !(state.norm_residual < norm_residual_prev) 
                converged_failed = true
                break
            end

            #Reset partstates
            norm_residual_prev = state.norm_residual
            state.partstates .= deepcopy(state0.partstates)
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

    determine_solvermode!(solver, state)

    #Recalculate f-start for this timestep
    fill!(system_arrays, 0.0)
    assemble_fstar!(globaldata.dh, state, globaldata)
    state.fˢ .= system_arrays.fⁱ

    return true
end

function newton_done(residual, Δg, ΔL, tol, solvermode) 
    if residual < tol 
        if solvermode == DISSIPATION
            if norm(Δg-ΔL) < tol
                return true
            else
                return false
            end
        end
        return true
    end
    return false
end

function _solve_dissipation_system(Kₜ, rₜ, q, fᵉ, fˢ, Δg, λ0, ΔL, solvermode)
    local ΔΔd, ΔΔλ
    if solvermode == DISSIPATION
        w = 0.0
        h = 0.5*(λ0*q - fˢ)
        #@show norm(Kₜ), norm(rₜ), norm(q), norm(fᵉ), norm(fˢ), norm(Δg), norm(ΔL), norm(λ0), solvermode
        KK = vcat(hcat(Kₜ, -q), hcat(h', w))
        ff = vcat(rₜ, -(Δg - ΔL))

        aa = KK\ff
        ΔΔd, ΔΔλ = (aa[1:end-1], aa[end])
    elseif solvermode == INCREMENT
        ΔΔd  = Kₜ\rₜ
        ΔΔλ = 0.0
    else
        error("No mode")
    end

    return ΔΔd, ΔΔλ
end

function set_initial_guess!(solver::DissipationSolver, state::StateVariables, ΔP, ⁿΔP, ntries)

    local ΔP_max, ΔP_min

    if state.solvermode == DISSIPATION
        ΔP_max = solver.ΔL_max
        ΔP_min = solver.ΔL_min
    elseif state.solvermode == INCREMENT
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
            ΔP *= (0.5^(0.25*(state.newton_itr-solver.optitr)))
        else
            #If the previous newton loop failed (ie ntries != 0),
            # half the step size
            ΔP /= 2
        end
    end

    if ΔP < ΔP_min
        ΔP = ΔP_min
    elseif ΔP > ΔP_max
        ΔP = ΔP_max
    end

    factor = ΔP / ⁿΔP

    if state.solvermode == DISSIPATION
        state.ΔL = ΔP
        state.Δλ *= factor
        state.Δd *= factor 
    elseif state.solvermode == INCREMENT
        state.Δλ = ΔP
        state.Δd *= factor
    end

    state.λ += state.Δλ
    state.d += state.Δd 
    
    return ΔP
end

function determine_solvermode!(solver::DissipationSolver, state::StateVariables)

    #Note: first step ΔL/abs(Δλ) == NaN, and will choose Incremenet
    #=if state.ΔL/abs(state.Δλ) < 1e-4
        return solver.solver_mode[] = INCREMENT
    else
        return solver.solver_mode[] = DISSIPATION
    end=#


    if state.solvermode == INCREMENT
        if state.ΔL >= solver.sw2d
            state.solvermode = DISSIPATION
            state.newton_itr = solver.optitr
        end
    elseif state.solvermode == DISSIPATION
        if state.ΔL/abs(state.Δλ) < solver.sw2i
            state.solvermode = INCREMENT
        end
    end


end