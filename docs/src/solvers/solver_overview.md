
# Solver overview

All solvers follow a similar pattern. Given a `state::StateVariables`, the solver tries to advance the 
state to the next timestep/loadstep. The step size is determined by the solver itself based on its
input parameters. If the current step diverges, most solvers tries to modify (e.g. decrease) the step size in order
to hopefully achieve convergence more easily.

The functions a solver (`<: AbstractSolver`) should implement are presented below. 

```@docs
Five.step!
Five.isdone
Five.should_abort
```