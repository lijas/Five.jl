export DissipationSolver

const ELASTIC = MODE1
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
    maxitr_first_step::Int = maxitr
    finish_criterion::Function = finish_criterion
end

function Base.isdone(solver::DissipationSolver, state::StateVariables, globaldata)
    return solver.finish_criterion(solver, state)
end

function should_abort(solver::DissipationSolver, ΔP, solvermode)
    if solvermode == ELASTIC && ΔP <= solver.Δλ_min 
        return true
    elseif solvermode == DISSIPATION && ΔP <= solver.ΔL_min 
        return true
    end
    return false
end

function step!(solver::DissipationSolver, state::StateVariables, globaldata)

    λ0     = state.λ
    q      = state.system_arrays.q
    state0 = deepcopy(state)

    #ΔP represents either ΔL or Δλ depending 
    # on the solver mode (INCREMENT or DISSIPATION).
    # Set ΔP to the privious time steps solution
    ΔP0 = ΔP = state.ΔL
    
    #
    if state.step == 1
        success = do_first_newton_step!(solver, state, globaldata)
        return success
    end

    converged_failed = true
    ntries = 0
    Δg = 0.0
    
    while converged_failed 
        ΔP = set_initial_guess!(solver, state, ΔP, ΔP0, ntries)
        
        println("Mode: $(state.solvermode), ntries: $(ntries), prev_λ = $(state0.λ-state0.Δλ), Δλ = $(state.Δλ), ΔL = $(state.ΔL)")
        
        state.newton_itr = 0
        ntries += 1
        while true
            state.newton_itr += 1
            fill!(state.system_arrays, 0.0)

            #Get internal force                                                                       
            assemble_stiffnessmatrix_and_forcevector!(globaldata.dh, state, globaldata)
            apply_constraints!(globaldata.dh, globaldata.constraints, state,  globaldata)
            
            #Normal stiffness matrix
            Kₜ = state.system_arrays.Kⁱ - state.system_arrays.Kᵉ
            rₜ = state.λ*q + state.system_arrays.fᵉ - state.system_arrays.fⁱ
            
            #Dissipaiton
            state.system_arrays.G[] = Δg = 1/2 * dot(state.Δd, state0.λ*q - state0.system_arrays.fᴬ)
            h = 0.5*(λ0*q - state0.system_arrays.fᴬ)
            w = 0.0
            if state.solvermode == ELASTIC
                state.system_arrays.G[] = 0.5*dot(state.Δλ*state0.d + state0.λ*state.Δd, q)
                w = 0.5*dot(state0.d,q)
                h = 0.5*state0.λ*q
            end
            
            apply_zero!(Kₜ, rₜ, globaldata.dbc)

            ΔΔd, ΔΔλ = _solve_dissipation_system(solver, Kₜ, rₜ, q, h, w, state.system_arrays.G[], state.ΔL)
            
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
            @printf("--->Newton %i, normr: %.5e, dissi: %.5e, Δλ: %.5e, Δg:%.5e \n", state.newton_itr, state.norm_residual, Δg, state.Δλ, state.system_arrays.G[])

            if newton_done(solver, state.norm_residual, state.system_arrays.G[], state.ΔL, state.solvermode)
                converged_failed = false
                break
            end

            maxitr = (state.step == 1) ? (solver.maxitr_first_step) : solver.maxitr
            if state.newton_itr >= maxitr || state.norm_residual > solver.max_residual || (abs(state.Δλ)/(solver.λ_max - solver.λ_min)) > 0.2 
                converged_failed = true
                break
            end


            #Reset partstates
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
    state.ΔL = state.system_arrays.G[]
    state.L += state.system_arrays.G[]

    determine_solvermode!(solver, state, Δg)

    #Recalculate f-start for this timestep
    assemble_fstar!(globaldata.dh, state, globaldata) #Stores in fᴬ
    #state.system_arrays.fᴬ .= state0.system_arrays.fᴬ

    return true
end

function newton_done(solver::DissipationSolver, residual, Δg, ΔL, solvermode) 
    if residual < solver.tol 
        if norm(Δg-ΔL) < solver.tol
            return true
        else
            return false
        end
    end
    return false
end

function _solve_dissipation_system(solver::DissipationSolver, Kₜ, rₜ, q, h, w, Δg, ΔL)
    local ΔΔd, ΔΔλ
    
    KK = vcat(hcat(Kₜ, -q), hcat(h', w))
    ff = vcat(rₜ, -(Δg - ΔL))

    aa = KK\ff
    ΔΔd, ΔΔλ = (aa[1:end-1], aa[end])

    return ΔΔd, ΔΔλ
end

function set_initial_guess!(solver::DissipationSolver, state::StateVariables, ΔP, ⁿΔP, ntries)

    local ΔP_max, ΔP_min

    if state.solvermode == DISSIPATION
        ΔP_max = solver.ΔL_max
        ΔP_min = solver.ΔL_min
    elseif state.solvermode == ELASTIC
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

function determine_solvermode!(solver::DissipationSolver, state::StateVariables, Δg)

    #Note: first step ΔL/abs(Δλ) == NaN, and will choose Incremenet
    #=if state.ΔL/abs(state.Δλ) < 1e-4
        return solver.solver_mode[] = INCREMENT
    else
        return solver.solver_mode[] = DISSIPATION
    end=#


    if state.solvermode == ELASTIC
        if Δg >= solver.sw2d
            state.solvermode = DISSIPATION
            state.newton_itr = solver.optitr
        end
    elseif state.solvermode == DISSIPATION
        if Δg/abs(state.Δλ) < solver.sw2i
            state.solvermode = ELASTIC_LOCAL
        end
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

        state.norm_residual = norm(r[JuAFEM.free_dofs(globaldata.dbc)])

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

    state.newton_itr = solver.optitr

    state.Δλ = solver.Δλ0
    state.λ = solver.Δλ0
    state.ΔL = dot(state.λ*state.Δd, state.system_arrays.q)
    state.L += state.ΔL

    @show state.ΔL

    return !conv_failed
end