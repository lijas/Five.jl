export ExplicitSolver, Output, solvethis

@with_kw struct ExplicitSolver{T} <: AbstractSolver{T} 
    Δt0::T
    Δt_min::Float64 = Δt0/1000
    Δt_max::Float64 = Δt0*1000
end


function Base.isdone(solver::ExplicitSolver, state::StateVariables, globaldata)
    return state.t >= globaldata.tend
end

function should_abort(solver::ExplicitSolver, state)
    return false
end

function step!(solver::ExplicitSolver, state, globaldata) 

    #
    M = state.system_arrays.Mᵈⁱᵃᵍ
    ⁿfⁱ = copy(state.system_arrays.fⁱ)
    ⁿfᵉ = copy(state.system_arrays.fᵉ)

    #Time update
    tᵢ = state.t
    tᵢ₊₀₅ = state.t + solver.Δt0/2
    tᵢ₊₁ = state.t + solver.Δt0
    state.t = tᵢ₊₁
    Δtᵢ₊₀₅  = solver.Δt0

    #Partial update of vel
    vᵢ₊₀₅ = state.v .+ (tᵢ₊₀₅  - tᵢ).*state.a 

    #Enforce veloctiy bc
    update!(globaldata.dbc, tᵢ₊₀₅ )
    apply!(vᵢ₊₀₅ , globaldata.dbc)

    #Update nodal disp
    state.Δd .= Δtᵢ₊₀₅ .* vᵢ₊₀₅ 
    state.d .+= state.Δd 

    @show tᵢ₊₁
    @show maximum(abs.(state.d))
    @show any(isnan.(state.d))

    fill!(state.system_arrays, 0.0)

    #Get internal force                                                                       
    @timeit "Assembling" assemble_forcevector!(globaldata.dh, state, globaldata)
    
    @timeit "ExternalForces" apply_external_forces!(globaldata.dh, globaldata.efh, state, globaldata)
    @timeit "Apply constraint" apply_constraints!(globaldata.dh, globaldata.constraints, state, globaldata)

    #Compute aᵢ₊₁ 
    update!(globaldata.dbc,tᵢ₊₁)
    
    #Assabmle mass matrix each iteration due to rigid bodies have non-constant mass matrix
    # @timeit "Update massmatrix" update_massmatrix!(dh, parts, materials, cellstates, elementinfos, M, dᵢ₊₁, vᵢ₊₀₅)
    f =state.system_arrays.fᵉ - state.system_arrays.fⁱ
    #apply_zero!(M, f, globaldata.dbc)
    @timeit "Solve" state.a .= M\f

    #Second partial update nodal vel
    state.v  = vᵢ₊₀₅  + (tᵢ₊₁ - tᵢ₊₀₅ )*state.a
    
    #Enforce veloctiy bc
    update!(globaldata.dbc, tᵢ₊₁ )
    apply!(state.v , globaldata.dbc)
    
    #Check energy balance 
    state.Wⁱ = 0.5*(state.Δd)'*(ⁿfⁱ + state.system_arrays.fⁱ)
    state.Wᵉ = 0.5*(state.Δd)'*(ⁿfᵉ + state.system_arrays.fᵉ)
    state.Wᵏ = 0.5*state.v'*M*state.v 

    return true

end