# # Bar example

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

material1 = LinearElastic(
    E = 210.0,
    ν = 0.3
)

bar1 = Part(
    material = material1,
    cellset = [1, 2],
    element = Five.BarElement{2}(
        area = 1.0,
    )
)
push!(data.parts, bar1)

midbar = Part(
    material = material1,
    cellset = [3],
    element = Five.BarElement{2}(
        area = 1.0/2,
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
    func = (X,t) -> 1.0
)
push!(data.external_forces, force)

solver = ArcLengthSolver(
    Δλ0 = -1.0,
    
    λ_max = 40.0,
    λ_min = -40.0,

    ΔL_max = 5.0,
    ΔL_min = 0.01,

    tol = 1e-4,
    maxsteps = 30,
    optitr = 10,
    maxitr = 20
)

state, data = build_problem(data)

result = solvethis(solver, state, data)

u = getproperty.(result.outputdata["reactionforce"].data, :displacement)
f = getproperty.(result.outputdata["reactionforce"].data, :fint)

using Test
@test isapprox( last(u)[1], -82.2, atol = 1e-1)
@test isapprox( last(f)[1], 9.3, atol = 1e-1)
# plot(u,f, mark=:o)
