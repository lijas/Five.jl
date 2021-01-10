export ArcLengthSolver

@with_kw struct ArcLengthSolver{T} <: AbstractSolver{T}
    Δλ0::T = 1.0
    λ_max::T
    λ_min::T

    ΔL0::T
    ΔL_min::T = ΔL0/10
    ΔL_max::T = ΔL0*10
    
    ψ::T = 0.0
    maxsteps::Int = 100

    #Newton
    tol::T = 1.0e-9
    max_residual::T = 1e12
    optitr::Int = 5
    maxitr::Int = 10
end

function Base.isdone(solver::ArcLengthSolver, state::StateVariables, globaldata)
    return state.λ > solver.λ_max || state.λ < solver.λ_min || state.step > solver.maxsteps
end

function should_abort(solver::ArcLengthSolver, ΔL, solvermode)
    return ΔL <= solver.ΔL_min 
end

function step!(solver::ArcLengthSolver, state::StateVariables, globaldata)
    λ0     = state.λ
    q      = state.system_arrays.q

    state0 = deepcopy(state)
    converged = false
    ntries = 0
    ΔL = (state.step == 1) ? solver.ΔL0 : state.ΔL

    δuₜ = zeros(ndofs(globaldata.dh))

    while !converged
        ΔL = set_initial_guess!(solver, state, globaldata, ΔL, ntries)
        
        state.newton_itr = 0
        ntries += 1
        
        println("Mode: ntries: $(ntries), λ = $(state.λ), Δλ = $(state.Δλ), ΔL = $(state.ΔL)")
        while true
            state.newton_itr += 1

            fill!(state.system_arrays, 0.0)

            #Get internal force                                                                       
            assemble_stiffnessmatrix_and_forcevector!(globaldata.dh, state, globaldata)
            apply_constraints!(globaldata.dh, globaldata.constraints, state,  globaldata)
            
            #Solve
            Kₜ = state.system_arrays.Kⁱ - state.system_arrays.Kᵉ
            rₜ = state.system_arrays.fⁱ - state.system_arrays.fᵉ - state.λ*q
            
            apply_zero!(Kₜ, rₜ, globaldata.dbc)
            
            state.detK = det(Kₜ)
            
            δuₜ .= (Kₜ\q)
            δū  = -(Kₜ\rₜ)

            #Solve arcleanth eq
            ψ = solver.ψ
            α1 = dot(δuₜ,δuₜ) + ψ*dot(q,q)
            α2 = 2*dot( (state.Δd + δū) , δuₜ) + 2*ψ^2 * state.Δλ * dot(q,q)
            α3 = dot( (state.Δd + δū), (state.Δd + δū) ) + ψ^2 * state.Δλ^2*dot(q,q) - ΔL^2

            _p = α2/α1; _q = α3/α1; 
            local δλ1, δλ2
            try 
                δλ1 = -_p/2 + sqrt((_p/2)^2 - _q)
                δλ2 = -_p/2 - sqrt((_p/2)^2 - _q)
            catch DomainError
                converged = false
                break
            end

            #Choose equation
            δu₁ = δū + δλ1*δuₜ
            δu₂ = δū + δλ2*δuₜ

            prev_Δd = state.step == 1 ? state.Δd : state0.Δd
            θ₁ = dot(state.Δd + δu₁, prev_Δd)
            θ₂ = dot(state.Δd + δu₂, prev_Δd) 

            δλ = θ₁ > θ₂ ? δλ1 : δλ2

            #Update d
            δd = δū + δλ*δuₜ
            state.Δd += δd
            state.Δλ += δλ
            state.d  += δd
            state.λ  += δλ

            #Check convergance
            state.norm_residual = norm(rₜ[JuAFEM.free_dofs(globaldata.dbc)])
            println("------>Normg: $(state.norm_residual), Δλ = $(state.Δλ)")
            
            #Check convergence
            if state.norm_residual < solver.tol
                converged = true
                break
            end

            if state.newton_itr >= solver.maxitr
                converged = false
                break
            end

            #Reset partstates
            state.partstates .= deepcopy(state0.partstates)
        end

        if !converged 
            if should_abort(solver, ΔL, state.solvermode)
                return false
            else
                #Reset the state
                copy!(state, state0)
            end
        end
    end
    
    state.system_arrays.fᴬ .= δuₜ #reuse fᴬ
    state.prev_detK = state0.detK

    return true
end

function set_initial_guess!(solver::ArcLengthSolver, state::StateVariables, globaldata, ΔL, ntries)

    prev_newton_itr = state.newton_itr
    detK = state.detK
    detK0 = state.prev_detK
    δuₜ = copy(state.system_arrays.fᴬ) #Reuse fs for ut
    Δλ0 = state.Δλ
    ⁿΔL = ΔL

    if state.step == 1
        prev_newton_itr = solver.optitr
        ⁿΔL = ΔL = solver.ΔL0
        Δλ0 = solver.Δλ0

        Kₜ = copy(state.system_arrays.Kⁱ - state.system_arrays.Kᵉ)
        _q = copy(state.system_arrays.q)
        
        apply_zero!(Kₜ, _q, globaldata.dbc)
        δuₜ .= Kₜ\_q

        detK = detK0 = det(Kₜ)
    end

    if ntries == 0 
        ΔL *= (0.5^(0.25*(prev_newton_itr-solver.optitr)))
    else
        ΔL /= 2
    end

    _sign = (sign(detK0) == sign(detK)) ? sign(Δλ0) : -sign(Δλ0)
    
    if ΔL < solver.ΔL_min
        ΔL = solver.ΔL_min
    elseif ΔL > solver.ΔL_max
        ΔL = solver.ΔL_max
    end

    state.ΔL = ΔL
    state.Δλ = _sign*ΔL/sqrt(dot(δuₜ,δuₜ))
    state.Δd = state.Δλ*δuₜ

    state.λ += state.Δλ
    state.d += state.Δd 
    
    return ΔL
end
