using Five

data = ProblemData(
    dim = 2,
    tend = 2.0
)

data.grid = generate_grid(QuadraticQuadrilateral, (5,5), Vec((0.0, 0.0)), Vec((1.0, 1.0)))

addvertexset!(data.grid, "botright", (x) -> x[1] == 10.0 && x[2] == 1.0)

material = MatLinearElastic(
    rho = 1.0,
    E = 1e5,
    nu = 0.3
)

#=@addmat MatLinearElastic, "Steel", 1 begin
    E = 1e5,
    nu = 0.3   
end=#

con1 = Dirichlet(
    set = getfaceset(data.grid, "top"),
    func = (x,t) -> (0.0, 1.0*t),
    field = :u,
    dofs = [1,2]
)
push!(data.dirichlet, con1)

part = Part{2,Float64}(
    element = SolidElement{2,2,RefCube,Float64}(
        celltype = JuAFEM.QuadraticQuadrilateral,
        qr_order = 4
    ),
    material = material,
    cellset = collect(1:getncells(data.grid))
)
push!(data.parts, part)

data.output[] = Output(
    interval = 0.1,
    runname = "contactexample",
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

contact = 
StaticPlaneContact{2}(
    slaves = getnodeset(data.grid, "bottomnodes")
    n = Vec{2}((0.0, 1.0)),
    a = Vec{2}((1.0, 0.0)),
    xm = Vec{2}((0.0, 0.0)),
    μ = 0.3,
    c_T = 1e5,
    penalty = 1e5
)
problem.contact[] = contact
    

solver = NewtonSolver(
    Δt0 = 0.1,
    Δt_max = 0.1,
)

state, data = build_problem(data)

solvethis(solver, state, data)
    
