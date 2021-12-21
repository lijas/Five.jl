
# Crisfield arc-length solver

This solver is based on the paper from 1981 by Crisfield (A fast incremental ...)

The arc-length constraint is defined as

```math
    \varphi(\boldsymbol a, \lambda) = \Delta \boldsymbol a \Delta \boldsymbol a^T + \Delta\lambda^2 - \Delta L^2 = 0
```

where $\Delta L$ is a user-defined parameter limiting the step size.