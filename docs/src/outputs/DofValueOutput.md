# DofValueOutputs

A `DofValueOutput` return the displacements and internal/external forces of the dofs in the specified set.

For example:

```@julia
    output = OutputData(
        type = DofValueOutput(
            field = :u,
            dofs = [2]
        ),
        interval = 0.1,
        set = getfaceset(data.grid, "bottom")
    )   
```
