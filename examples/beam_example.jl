using Five

data = ProblemData(
    dim = 2,
    tend = 1.0
)

data.grid = generate_grid(Quadrilateral, (10,1), Vec((0.0, 0.0)), Vec((10.0, 1.0)))

addvertexset!(data.grid, "topright", (x) -> x[1] == 10.0 && x[2] == 1.0)

material = MatNeoHook(
    E = 1e5,
    ν = 0.3
)

#=@addmat MatLinearElastic, "Steel", 1 begin
    E = 1e5,
    nu = 0.3   
end=#

con1 = Dirichlet(
    set = getfaceset(data.grid, "left"),
    func = (x,t) -> (0.0, 0.0),
    field = :u,
    dofs = [1,2]
)
push!(data.dirichlet, con1)

part = Part{2,Float64}(
    element = Five.SolidElementQuad(),
    material = material,
    cellset = collect(1:getncells(data.grid))
)
push!(data.parts, part)

data.output[] = Output(
    interval = 0.1,
    runname = "Beamexample",
    savepath = "."
)

output = OutputData(
    type = DofValueOutput(
        field = :u,
        dofs = 1:2
    ),
    interval = 0.1,
    set = getfaceset(data.grid, "left")
)
data.outputdata["reactionforce"] = output

vtkoutput = VTKCellOutput(
    type = MaterialStateOutput(
        field = :σ
    ),
    func = mean,
)
Five.push_vtkoutput!(data.output[], vtkoutput)

force = PointForce(
    field = :u,
    comps = [2],
    set = getvertexset(data.grid, "topright"),
    func = (X,t) -> -10.0*t
)
push!(data.external_forces, force)

solver = NewtonSolver(
    Δt0 = 0.1,
    Δt_max = 0.1,
)

state, data = build_problem(data)

solvethis(solver, state, data)
    
