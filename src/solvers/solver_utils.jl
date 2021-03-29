

#Some common stuff for arc-length solvers
finish_criterion(solver::AbstractSolver, state) = state.λ > solver.λ_max || state.λ < solver.λ_min || state.step > solver.maxsteps
finish_criterion1(solver::AbstractSolver, state::StateVariables) = (state.Δλ < 0 && state.λ < solver.λ_min) || state.step > solver.maxsteps
finish_criterion2(solver::AbstractSolver, state::StateVariables) = (state.Δλ > 0 && state.λ > solver.λ_max) || state.step > solver.maxsteps
