export ArcLengthSolver

@with_kw struct ArcLengthSolver{T} <: AbstractSolver{T}
    Δλ0::T = 1.0
    λ_max::T
    λ_min::T

    ΔL_min::T
    ΔL_max::T
    
    ψ::T = 0.0
    maxsteps::Int = 100

    #Newton
    tol::T = 1.0e-9
    max_residual::T = 1e12
    optitr::Int = 5
    maxitr::Int = 10
    maxitr_first_step::Int = maxitr
    finish_criterion::Function = finish_criterion
end


"""
    ArcLengthSolver{T}

Arguments:
    - `Δλ0`: Initial size of load parameter
    - `λ_max`: 
    - `λ_min`: 
    - `ψ`: 
    - `maxsteps`: Maximum number of steps
    - `max_residual`: Maximum value of the residual before aborting a step
    - `optitr`: 
    - `maxitr`: 
    - `maxitr_first_step`: 
    - `finish_criterion`: Function 
"""
ArcLengthSolver

function Base.isdone(solver::ArcLengthSolver, state::StateVariables, globaldata)
    return solver.finish_criterion(solver, state)
end

function should_abort(solver::ArcLengthSolver, state::StateVariables, globaldata)
    return state.ΔL <= solver.ΔL_min 
end

function step!(solver::ArcLengthSolver, state::StateVariables, globaldata, ntries::Int = 0)
    ch = globaldata.dbc
    dh = globaldata.dh
    partstates0 = deepcopy(state.partstates)
    Δd0 = copy(state.v) #Use v (velocity) as displacement increment

    q  = state.system_arrays.q
    δuₜ = copy(state.system_arrays.fᴬ) #Realias fᴬ

    if state.step == 1
        Kₜ = state.system_arrays.Kⁱ
        _q = copy(q)
        
        apply!(Kₜ, _q, ch, true; strategy = Ferrite.APPLY_TRANSPOSE)
        δuₜ .= Kₜ\_q
        apply_zero!(δuₜ, ch)

        state.ΔL = abs(solver.Δλ0) * sqrt(dot(δuₜ,δuₜ))
        state.newton_itr = solver.optitr
        state.Δλ = solver.Δλ0
        state.system_arrays.fᴬ .= δuₜ
        state.detK  = det(Kₜ)
    end

    set_initial_guess!(solver, state, ntries)

    ΔL = state.ΔL
    detK0 = state.detK
    Δd = state.v #Realias the velocity to Δd
    println("Step: $(state.step), λ: $(state.λ), Δλ: $(state.Δλ), ΔL: $(state.ΔL)")

    state.newton_itr = 0
    while true
        state.newton_itr += 1

       zero_out_systemarrays!(state.system_arrays)

        @timeit "Assembling"       assemble_stiffnessmatrix_and_forcevector!(dh, state, globaldata)
        @timeit "Apply constraint" apply_constraints!(dh, globaldata.constraints, state, globaldata)
        
        #Solve
        Kₜ = state.system_arrays.Kⁱ - state.system_arrays.Kᵉ
        rₜ = state.system_arrays.fⁱ - state.system_arrays.fᵉ - state.λ*q

        apply_zero!(Kₜ, rₜ, globaldata.dbc)

        δuₜ .= (Kₜ\q)
        δū  = -(Kₜ\rₜ)

        #Solve arcleanth eq
        ψ = solver.ψ
        α1 =  dot(δuₜ,δuₜ)  + ψ*dot(q,q)
        α2 =2*dot( (Δd + δū) , δuₜ)  + 2*ψ^2 * state.Δλ * dot(q,q)
        α3 =  dot( (Δd + δū), (Δd + δū) )  + ψ^2 * state.Δλ^2*dot(q,q) - ΔL^2

        _p = α2/α1; _q = α3/α1; 
        local δλ1, δλ2
        try 
            δλ1 = -_p/2 + sqrt((_p/2)^2 - _q)
            δλ2 = -_p/2 - sqrt((_p/2)^2 - _q)
        catch DomainError
            return false
        end

        #Choose equation
        δu₁ = δū .+ δλ1*δuₜ
        δu₂ = δū .+ δλ2*δuₜ

        θ₁ = dot(Δd + δu₁, state.step == 1 ? Δd : Δd0)
        θ₂ = dot(Δd + δu₂, state.step == 1 ? Δd : Δd0) 
        #@show (δλ1, θ₁), (δλ2, θ₂)
        
        δλ = θ₁ > θ₂ ? δλ1 : δλ2

        #Update d
        δd = δū .+ δλ*δuₜ
        state.d  .+= δd
        state.Δλ += δλ
        Δd .+= δd
        state.λ  += δλ

        #Check convergance
        if norm(state.λ*q) <= 1e-10
            state.norm_residual = norm(rₜ[Ferrite.free_dofs(globaldata.dbc)])
        else
            state.norm_residual = norm(rₜ[Ferrite.free_dofs(globaldata.dbc)])/norm(state.λ*q)
        end
        println("------>Normg: $(state.norm_residual), λ: $(state.λ), Δλ = $(state.Δλ)")

        maxitr = (state.step == 1) ? (solver.maxitr_first_step) : solver.maxitr
        if state.newton_itr >= maxitr || state.norm_residual > solver.max_residual
            return false
        end

        if state.norm_residual < solver.tol
            break
        end

        state.partstates .= deepcopy(partstates0)
    end

    state.system_arrays.fᴬ .= δuₜ #reuse fᴬ

    return true
end

function set_initial_guess!(solver::ArcLengthSolver, state::StateVariables, ntries)
 
    δuₜ = copy(state.system_arrays.fᴬ)
    detK0 = state.detK
    state.detK = det(state.system_arrays.Kⁱ) #This assumes that the stiffness matrix is kept untouched between steps.

    if ntries == 0 
        state.ΔL *= (0.5^(0.25*(state.newton_itr - solver.optitr)))
    else
        state.ΔL /= 2^(ntries-1)
    end

    #Check if determinant has changed sign
    _sign = (sign(state.detK) == sign(detK0)) ? sign(state.Δλ) : -sign(state.Δλ)

    if state.ΔL < solver.ΔL_min
        state.ΔL = solver.ΔL_min
    elseif state.ΔL > solver.ΔL_max
        state.ΔL = solver.ΔL_max
    end

    state.Δλ = _sign*state.ΔL/sqrt(dot(δuₜ,δuₜ))
    state.v .= state.Δλ*δuₜ # ./Δt

    state.λ += state.Δλ
    state.d .+= state.v # *.Δt 
    
end
