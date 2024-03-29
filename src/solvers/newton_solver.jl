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

    linearsolver::Optional{LinearSolve.SciMLLinearSolveAlgorithm} = nothing
    preconditioner = IdentityProconditioner
end

function Base.isdone(solver::NewtonSolver, state::StateVariables, globaldata)
    return state.t >= globaldata.tend
end

function should_abort(solver::NewtonSolver, state, globaldata)
    return state.Δt <= solver.Δt_min
end

function step!(solver::NewtonSolver, state::StateVariables, globaldata, ntries=0)

    ch = globaldata.dbc
    dh = globaldata.dh
    partstates0 = deepcopy(state.partstates)
    d0 = copy(state.d)

    set_initial_guess(solver, state, ntries)

    #Apply boundary conditions
    update!(ch, state.t)
    apply!(state.d, ch)
    state.newton_itr = 0
    while true
        state.newton_itr += 1

        fill!(state.system_arrays, 0.0)

        @timeit "Assembling"       assemble_stiffnessmatrix_and_forcevector!(dh, state, globaldata)
        @timeit "ExternalForces"   apply_external_forces!(dh, globaldata.efh, state, globaldata)
        @timeit "Apply constraint" apply_constraints!(dh, globaldata.constraints, state, globaldata)

        r = state.system_arrays.fⁱ - state.system_arrays.fᵉ
        K = state.system_arrays.Kⁱ #- state.system_arrays.Kᵉ
        #Solve 
        apply_zero!(K, r, ch)
        
        #Kt = ThreadedSparseMatrixCSC(K)' #Note the transpose
        prob = LinearSolve.LinearProblem(K, -r) 
        precon = solver.preconditioner(K)
        linsolve = LinearSolve.init(prob, solver.linearsolver, Pl=precon)

        @timeit "Solve" sol = LinearSolve.solve(linsolve)
        ΔΔd = sol.u
        apply_zero!(ΔΔd, ch)
        
        state.norm_residual = norm(r)

        state.d .+= ΔΔd
        state.v .+= ΔΔd # (state.d .- d0)#./state.Δt

        @info "[NEWTONSOLVER] ----resdiual: $(state.norm_residual)"

        maxitr = (state.step == 1) ? (solver.maxitr_first_step) : solver.maxitr
        if state.newton_itr >= maxitr || state.norm_residual > solver.max_residual
            return false
        end

        if state.norm_residual < solver.tol
            break
        end

        state.partstates .= deepcopy(partstates0)
    end
    return true
end

function set_initial_guess(solver::NewtonSolver,  state::StateVariables, ntries::Int)

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
    state.d += state.v#(state.v*Δt) * factor
end