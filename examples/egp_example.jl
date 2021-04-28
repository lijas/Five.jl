using Five

tmax = 500.0
data = ProblemData(
    dim = 3,
    tend = tmax
)

data.grid = generate_grid(Hexahedron, (1,1,1), Vec((0.0, 0.0, 0.0)), Vec((1.0, 1.0, 1.0)))

addvertexset!(data.grid, "leftvert", (x) -> x[1] == 0.0 && x[2] == 0.0 && x[3] == 0.0)

material = MatEGP()

con1 = Dirichlet(
    set = getfaceset(data.grid, "left"),
    func = (x,t) -> (0.0),
    field = :u,
    dofs = [1]
)
push!(data.dirichlet, con1)

con1 = Dirichlet(
    set = getvertexset(data.grid, "leftvert"),
    func = (x,t) -> (0.0,0.0),
    field = :u,
    dofs = [2,3]
)
push!(data.dirichlet, con1)

loadf(t) = begin
    if t < tmax/2
        return (-t/tmax) * 0.3
    else
        return (t-tmax)/tmax * 0.3
    end
end

con1 = Dirichlet(
    set = getfaceset(data.grid, "right"),
    func = (x,t) -> (loadf(t), ),
    field = :u,
    dofs = [1]
)
push!(data.dirichlet, con1)

part = Part{3,Float64}(
    element = SolidElement{3,1,RefCube,Float64}(
        celltype = JuAFEM.Hexahedron,
        qr_order = 3
    ),
    material = material,
    cellset = collect(1:getncells(data.grid))
)
push!(data.parts, part)

data.output[] = Output(
    interval = 100.0,
    runname = "egptest",
    savepath = "."
)

output = OutputData(
    type = DofValueOutput(
        field = :u,
        dofs = [1]
    ),
    interval = 1.0,
    set = getfaceset(data.grid, "right")
)
data.outputdata["reactionforce"] = output

solver = NewtonSolver(
    Δt0 = 0.1,
    Δt_max = 10.0,
    tol = 1e-4
)

state, data = build_problem(data)

output = solvethis(solver, state, data)
  
d = [output.outputdata["reactionforce"].data[i].displacement for i in 1:length(output.outputdata["reactionforce"].data)]
f = [output.outputdata["reactionforce"].data[i].fint for i in 1:length(output.outputdata["reactionforce"].data)]

using Plots; plotly()
plot(d,f)