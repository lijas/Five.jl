using Five

const α = pi/3;
const Δ = 200.0
const P = -30

data = ProblemData(
    dim = 2,
    tend = 1.0
)

function generate_bars()
    L = 0.5Δ/cos(α)
    nodecoords = [Vec(0.0,0.0), Vec(0.5Δ, 0.5Δ/tan(α)), Vec(Δ,0.0), Vec(0.5Δ, 0.5Δ/tan(α) + L)]

    nodes = [Node{2,Float64}(x) for x in nodecoords]
    cells = [Line2D((1,2)), Line2D((2,3)), Line2D((2,4))]
    grid = Grid(cells,nodes)
    
    addvertexset!(grid, "left", (x)-> x[1] ≈ 0.0)
    addvertexset!(grid, "right", (x)-> x[1] ≈ Δ)
    addvertexset!(grid, "topmid", (x)-> x[1] ≈ Δ/2 && x[2] ≈ 0.5Δ/tan(α)+L)
    addvertexset!(grid, "midmid", (x)-> x[1] ≈ Δ/2)# && x[2] ≈ 0.5Δ/tan(α)+L)
end

data.grid = generate_bars()

material1 = MatLinearElastic(
    E = 210.0,
    nu = 0.3
)

bar1 = Part{2,Float64}(
    material = material1,
    cellset = [1, 2],
    element = BarElement{2}(
        area = 1.0,
    )
)
push!(data.parts, bar1)

midbar = Part{2,Float64}(
    material = material1,
    cellset = [3],
    element = BarElement{2}(
        area = 1.0,
    )
)
push!(data.parts, midbar)

con = Dirichlet(
    set = getvertexset(data.grid, "left"),
    func = (x,t) -> (0.0, 0.0),
    field = :u,
    dofs = [1,2]
)
push!(data.dirichlet, con)

con = Dirichlet(
    set =  getvertexset(data.grid, "right") ,
    func = (x,t) -> (0.0, 0.0),
    field = :u,
    dofs = [1,2]
)
push!(data.dirichlet, con)

con3 = Dirichlet(
    set =  Set([first(getvertexset(data.grid, "midmid"))]),
    func = (x,t) -> (0.0),
    field = :u,
    dofs = [1,]
)
push!(data.dirichlet, con3)

con3 = Dirichlet(
    set =  Set([first(getvertexset(data.grid, "topmid"))]),
    func = (x,t) -> (0.0),
    field = :u,
    dofs = [1,]
)
push!(data.dirichlet, con3)


data.output[] = Output(
    runname = "barexample",
    savepath = ".",
    interval = 0.0,
)

output = OutputData(
    type = DofValueOutput(
        field = :u,
        dofs = [2,]
    ),
    interval = 0.00,
    set = getvertexset(data.grid, "topmid")
)
data.outputdata["reactionforce"] = output

#=output = VTKPointData(
    type = MaterialOutput(
        field = :σ,
    ),
    part = bars
)=#

force = PointForce(
    field = :u,
    comps = [2,],
    set = getvertexset(data.grid, "topmid"),
    func = (X,t) -> -1.0
)
push!(data.external_forces, force)

#=solver = NewtonSolver(
    Δt0 = 0.01,
    Δt_max = 0.01,
    Δt_min = 0.001,
    tol = 1e-4
)=#

output = OutputData(
    type = EigenOutput(),
    interval = 0.0,
    set      = Set()
)
data.outputdata["eigen"] = output

solver = ArcLengthSolver(
    Δλ0 = 1.0,
    λ_max = 40.0,
    λ_min = -40.0,
    ΔL_min = 0.01,
    ΔL_max = 5.0,
    tol = 1e-4,
    maxsteps = 60,
    optitr = 5
)

state, globaldata = build_problem(data)

output = solvethis(solver, state, globaldata)

u = getproperty.(output.outputdata["reactionforce"].data, :displacement)
f = getproperty.(output.outputdata["reactionforce"].data, :fint)
plot(u,f, mark=:o)



for imode in 1:7
    d = output.outputdata["eigen"].data[end].d
    vec = output.outputdata["eigen"].data[end].evec[:,imode]

    state.d .= vec*1.0
    Five._vtk_add_state!(output, state, globaldata, outputname = "eigenval$imode");
end

#=
for imode in 1:6
    vals = []
    for step in 1:length(output.outputdata["eigen"].data)
        push!(vals, output.outputdata["eigen"].data[step].eval[imode])
    end

    fig = plot(vals./norm(vals), mark = :o)
    plot!(fig, f./norm(f))
    display(fig)
end
=#
