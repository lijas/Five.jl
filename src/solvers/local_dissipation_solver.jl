export LocalDissipationSolver

#const INCREMENT_LOCAL = MODE1
#const DISSIPATION_LOCAL = MODE2

@with_kw struct LocalDissipationSolver{T} <: AbstractSolver{T}
    Δλ0::T
    Δλ_min::T = Δλ0/1000
    Δλ_max::T = Δλ0*1000
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
    maxitr_first_step::Int = maxitr
    finish_criterion::Function = finish_criterion
end

#Type
Base.@kwdef struct LocalDissipationSolverModes
    INCREMENT::SolverMode = MODE1 #Used first time step
    RIKS::SolverMode = MODE2
    DISSIPATION::SolverMode = MODE3
end
DissipationSolverModes = LocalDissipationSolverModes()

function Base.isdone(solver::LocalDissipationSolver, state::StateVariables, globaldata)
    return solver.finish_criterion(solver, state)
end

function should_abort(solver::LocalDissipationSolver, state::StateVariables, globaldata)
    if state.solvermode == DissipationSolverModes.RIKS
        return state.Δλ <= solver.Δλ_min    
    elseif state.solvermode == DissipationSolverModes.DISSIPATION
        return state.ΔL <= solver.ΔL_min    
    else
        error("Unreachable code")
    end
end

function step!(solver::LocalDissipationSolver, state::StateVariables, globaldata, ntries::Int = 0)

    λ0     = state.λ
    q      = state.system_arrays.q
    fᵉ     = state.system_arrays.fᵉ
    partstates0 = deepcopy(state.partstates)
    
    Δg = 0.0

    set_initial_guess!(solver, state, ntries)
    Δd0 = copy(state.v) #Needed for riks solver

    println("Mode: $(state.solvermode), ntries: $(ntries), prev_state.λ = $(state.λ-state.Δλ), Δλ = $(state.Δλ), ΔL = $(state.ΔL)")
    
    state.newton_itr = 0
    while true
        state.newton_itr += 1

        fill!(state.system_arrays, 0.0)
        
        @timeit "Calculate dissipation" assemble_dissipation!(globaldata.dh, state, globaldata)
        Δg = state.system_arrays.G[]
        
        #Get internal force                                                                       
        @timeit "Assembling" assemble_stiffnessmatrix_and_forcevector!(globaldata.dh, state, globaldata)
        @timeit "Apply constraint" apply_constraints!(globaldata.dh, globaldata.constraints, state, globaldata)
        
        #Normal stiffness matrix
        Kₜ = state.system_arrays.Kⁱ - state.system_arrays.Kᵉ
        rₜ = state.λ*q + state.system_arrays.fᵉ - state.system_arrays.fⁱ
        
        apply_zero!(Kₜ, rₜ, globaldata.dbc)

        @timeit "Solve system" ΔΔd, ΔΔλ, _success = _solve_dissipation_system(solver, Kₜ, rₜ, q, state.system_arrays.fᵉ, state.system_arrays.fᴬ, Δg, Δd0, λ0, state.ΔL, state.solvermode)
        
        if !_success
            @info "Failed here"
            return false
        end

        state.v += ΔΔd
        state.Δλ += ΔΔλ
        state.d  += ΔΔd
        state.λ  += ΔΔλ

        #Check convergance
        scaledtol  = solver.tol * max(norm(fᵉ), norm(state.λ*q), 1.0)
        state.norm_residual = norm(rₜ[Ferrite.free_dofs(globaldata.dbc)])

        println("------>Newton $(state.newton_itr): 
                 $(rpad("Residual: $(state.norm_residual),", 20)) 
                 $(rpad("Δg=$(Δg),", 20)) 
                 $(rpad("Δλ=$(state.Δλ),", 20)) 
                 $(rpad("λ=$(state.λ),", 20))")

        if newton_done(state.norm_residual, Δg, state.ΔL, state.solvermode, scaledtol, solver.tol)
            @info "convergence found"
            break
        end

        maxitr = (state.step == 1) ? (solver.maxitr_first_step) : solver.maxitr
        if state.newton_itr >= maxitr || state.norm_residual > solver.max_residual || (abs(state.Δλ)/(solver.λ_max - solver.λ_min)) > 0.2 
            @warn "Newton iterations failed"
            return false
        end

        #Reset partstates
        state.partstates .= deepcopy(partstates0)
    end


    #Update the dissipation increaments for this timestep
    state.ΔL = Δg
    state.L += Δg
    state.Δt = Δg
    state.t += Δg

    determine_solvermode!(solver, state)

    return true
end

function newton_done(residual, Δg, ΔL, solvermode, scaledtol::Float64, tol_dissipation::Float64) 
    if residual < scaledtol 
        if solvermode == DissipationSolverModes.DISSIPATION
            if norm(Δg-ΔL) < tol_dissipation
                return true
            else
                return false
            end
        end
        return true
    end
    return false
end

function _solve_dissipation_system(solver::LocalDissipationSolver, Kₜ, rₜ, q, fᵉ, fᴬ, Δg, Δd0, λ0, ΔL, solvermode)
    local ΔΔd, ΔΔλ
    if solvermode == DissipationSolverModes.INCREMENT
        ΔΔd  = Kₜ\rₜ
        ΔΔλ = 0.0
    elseif solvermode == DissipationSolverModes.DISSIPATION
        w = 0.0
        h = fᴬ
        φ = Δg - ΔL

        solvefull = 2
        if solvefull == 1
            KK = vcat(hcat(Kₜ, -q), hcat(h', w))
            ff = vcat(rₜ, -(Δg - ΔL))

            local aa
            try
                aa = KK\ff
            catch error
                error("Matrix not invertible")
                return 0.0 .* copy(h), 0.0, false
            end
            ΔΔd, ΔΔλ = (aa[1:end-1], aa[end])
        else
            RHS = hcat(rₜ, -q)
            SOL = Kₜ\RHS

            dᴵ  = SOL[:,1]
            dᴵᴵ = SOL[:,2]

            denom = (dot(h,dᴵᴵ) - w)
            ΔΔd = dᴵ - 1/(denom) * (dot(h,dᴵ) + φ)*dᴵᴵ
            ΔΔλ = -φ - 1/(denom) * (-dot(h,dᴵ) - φ*(1 + dot(h,dᴵᴵ) - w))

        end
    elseif solvermode == DissipationSolverModes.RIKS

        RHS = hcat(rₜ, -q)
        SOL = Kₜ\RHS

        dᴵ  = SOL[:,1]
        dᴵᴵ = SOL[:,2]     

        ΔΔλ = dot(Δd0,dᴵ)/dot(Δd0,dᴵᴵ)
        ΔΔd = -ΔΔλ*dᴵᴵ + dᴵ
    else
        @show solvermode
        error("No mode")
    end

    return ΔΔd, ΔΔλ, true
end

function set_initial_guess!(solver::LocalDissipationSolver, state::StateVariables, ntries)

    ΔP_max = solver.ΔL_max
    ΔP_min = solver.ΔL_min
    ⁿΔP = ΔP = state.ΔL
    if state.solvermode == DissipationSolverModes.INCREMENT || state.solvermode == DissipationSolverModes.RIKS
        ΔP_max = solver.Δλ_max
        ΔP_min = solver.Δλ_min
        ⁿΔP = ΔP = state.Δλ
    end

    if state.step == 1
        ⁿΔP = ΔP = solver.Δλ0 * (1/2)^ntries
    else
        if ntries == 0 
            #For the first try, increase/decrease the step 
            # based on the number of newton_iteration in the previeus 
            # converged solution
            ΔP *= (0.5^(0.1*(state.newton_itr-solver.optitr)))
            @show ΔP
        else
            #If the previous newton loop failed (ie ntries != 0),
            # half the step size
            ΔP *= (1/2)^ntries
        end
    end

    if ΔP < ΔP_min
        ΔP = ΔP_min
    elseif ΔP > ΔP_max
        ΔP = ΔP_max
    end

    factor = ΔP / ⁿΔP
    @show factor

    if state.solvermode == DissipationSolverModes.DISSIPATION
        state.ΔL = ΔP
        state.Δλ *= factor
        state.v *= factor 
    elseif state.solvermode == DissipationSolverModes.INCREMENT
        state.Δλ = ΔP
        state.v *= factor
    elseif state.solvermode == DissipationSolverModes.RIKS
        state.Δλ = ΔP
        state.v *= factor
    end

    state.λ += state.Δλ
    state.d .+= state.v 
end

function determine_solvermode!(solver::LocalDissipationSolver, state::StateVariables)

    #Note: first step ΔL/abs(Δλ) == NaN, and will choose Incremenet
    #=if state.ΔL/abs(state.Δλ) < 1e-4
        return solver.solver_mode[] = INCREMENT_LOCAL
    else
        return solver.solver_mode[] = DISSIPATION_LOCAL
    end=#

    #Only use the incremental solver in the first step
    if state.step == 1
        state.solvermode = DissipationSolverModes.RIKS
        state.newton_itr = solver.optitr #Dont assume anything about the number of newton iterations when swithing mode
    end

    if state.solvermode == DissipationSolverModes.RIKS
        if state.ΔL >= solver.sw2d
            state.solvermode = DissipationSolverModes.DISSIPATION
            state.newton_itr = solver.optitr
        end
    elseif state.solvermode == DissipationSolverModes.DISSIPATION
        if state.ΔL/abs(state.Δλ) < solver.sw2i
            state.solvermode = DissipationSolverModes.RIKS
        end
    else
        error("Unreachable")
    end


end