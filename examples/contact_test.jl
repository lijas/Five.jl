using Five

data = ProblemData(
    dim = 2,
    tend = 0.5
)

data.grid = generate_grid(Quadrilateral, (5,5), Vec((0.0, 0.0)), Vec((1.0, 1.0)))

addnodeset!(data.grid, "bottomnodes", (x) -> x[2] ≈ 0.0)

material = MatLinearElastic(
    rho = 1.0,
    E = 1e5,
    nu = 0.3
)

#=@addmat MatLinearElastic, "Steel", 1 begin
    E = 1e5,
    nu = 0.3   
end=#

function bcmove(x,t)

    if t < 0.1
        return (0.0, -1.0*t)
    else
        return (-1.0 * (t-0.1), -1.0*0.1)
    end

end

con1 = Dirichlet(
    set = getfaceset(data.grid, "top"),
    func = (x,t) -> bcmove(x,t),
    field = :u,
    dofs = [1,2]
)
push!(data.dirichlet, con1)

part = Part{2,Float64}(
    element = SolidElement{2,1,RefCube,Float64}(
        celltype = JuAFEM.Quadrilateral,
        qr_order = 4
    ),
    material = material,
    cellset = collect(1:getncells(data.grid))
)
push!(data.parts, part)

data.output[] = Output(
    interval = 0.0,
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
Five.StaticPlaneContact{2}(
    slaveset = getnodeset(data.grid, "bottomnodes") |> (x)->nodeset_to_vertexset(data.grid, x) |> collect,
    n = Vec{2}((0.0, 1.0)),
    a = Vec{2}((1.0, 0.0)),
    xm = Vec{2}((0.0, -0.03)),
    μ = 0.01,
    c_T = 1e4,
    penalty = 1e5
)
data.contact[] = contact
    

solver = NewtonSolver(
    Δt0 = 0.01,
    Δt_max = 0.01,
    Δt_min = 0.01,
)

state, data = build_problem(data)

solvethis(solver, state, data)
    
