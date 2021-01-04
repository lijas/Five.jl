using Five

const α = pi/4;
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
    cellset = [1],
    element = BarElement{2}(
        area = 1.0,
    )
)
push!(data.parts, bar1)

bar2 = Part{2,Float64}(
    material = material1,
    cellset = [2],
    element = BarElement{2}(
        area = 1.0,
    )
)
push!(data.parts, bar2)

midbar = Part{2,Float64}(
    material = material1,
    cellset = [3],
    element = BarElement{2}(
        area = 1.0,
    )
)
push!(data.parts, midbar)

con1 = Dirichlet(
    set = getvertexset(data.grid, "left"),
    func = (x,t) -> (0.0, 0.0),
    field = :u,
    dofs = [1,2]
)
push!(data.dirichlet, con1)

con1 = Dirichlet(
    set =  getvertexset(data.grid, "right") ,
    func = (x,t) -> (0.0, 0.0),
    field = :u,
    dofs = [1,2]
)
push!(data.dirichlet, con1)

con = Dirichlet(
    set = Set([first(getvertexset(data.grid, "midmid"))]),
    func = (x,t) -> (0.0),
    field = :u,
    dofs = [1,]
)
push!(data.dirichlet, con)

con2 = Dirichlet(
    set = Set([first(getvertexset(data.grid, "topmid"))]),
    func = (x,t) -> (0.0),
    field = :u,
    dofs = [1,]
)
push!(data.dirichlet, con2)

con3 = Dirichlet(
    set =  Set([first(getvertexset(data.grid, "topmid"))]),
    func = (x,t) -> (-100.0*t),
    field = :u,
    dofs = [2,]
)
#push!(data.dirichlet, con3)


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
    tol = 1e-9
)=#

solver = ArcLengthSolver(
    λ_max = 70.0,
    λ_min = -70.0,
    ΔL0 = 1.0,
    tol = 1e-8,
    maxsteps = 150
)

state, data = build_problem(data)

result = solvethis(solver, state, data)

u = getproperty.(result.outputdata["reactionforce"].data, :displacement)
f = getproperty.(result.outputdata["reactionforce"].data, :fint)


#=
function setup_test_example2()

    dim = 2
    T = Float64

    #Grid
    α = pi/3;
    Δ = 200.0
    P = -30

    nodecoords = [Vec(0.0,0.0), Vec(0.5Δ, 0.5Δ/tan(α)), Vec(Δ,0.0)]

    nodes = [Node{2,Float64}(x) for x in nodecoords]
    cells = [Cell{2,2,1}((1,2)), Cell{2,2,1}((2,3))]

    grid = Grid(cells,nodes)

    #
    parts = Five.AbstractPart{2}[]
    dbc   = JuAFEM.Dirichlet[]
    exfor = Five.AbstractExternalForce[]
    output = Dict{String, Five.AbstractOutput}()
    cnstr = Five.AbstractExternalForce[]

    #Part
    barinfo1 = Five.BarInfo(1.0)
    barinfo2 = Five.BarInfo(1.0)
    material = Five.MatLinearElastic{1}(0.0, 210.0, 0.3)
    part1 = Part{dim,T}(material, barinfo1, [1], Five.BarElement{dim}()) 
    part2 = Part{dim,T}(material, barinfo2, [2], Five.BarElement{dim}()) 

    push!(parts, part1)
    push!(parts, part2)

    #
    addvertexset!(grid, "left", (x)-> x[1] ≈ 0.0)
    addvertexset!(grid, "right", (x)-> x[1] ≈ Δ)
    addvertexset!(grid, "mid", (x)-> x[1] ≈ Δ/2)

    #
    dbc1 = JuAFEM.Dirichlet(:u, getvertexset(grid, "left"), (x,t)->[0.0, 0.0], [1,2])
    push!(dbc, dbc1)

    dbc1 = JuAFEM.Dirichlet(:u, getvertexset(grid, "right"), (x,t)->[0.0, 0.0], [1,2])
    push!(dbc, dbc1)

    dbc1 = JuAFEM.Dirichlet(:u, getvertexset(grid, "mid"), (x,t)->[t*-100.0], [2])
    #push!(dbc, dbc1)

    #
    force = VectorForce(BarElement{2}(), first(getvertexset(grid, "mid")), :u, [2], (t)->[-1.0])
    push!(exfor, force)

    #
    #force = VectorForce(BarElement{2}(), first(getvertexset(grid, "mid")), :u, [1], (t)->[0.0])
    #push!(cnstr, force)

    #
    rundir = "examples/paraview/pv1/"
    filename = "adap"

    #
    output["force"] = DofOutput(BarElement{2}(), first(getvertexset(grid, "mid")), :u,  2, 0.0)

    #Simulation time
    starttime = 0.0
    endtime = 1.0;
    dt = 0.01;
    plotinterval = 0.01#0.01;

    data = Five._build_problem(rundir, filename, grid, parts, dbc, exfor, output, cnstr)

    #return StaticSolver2{dim,T}(data, starttime, endtime, dt);
    
    return ArcLengthSolver{dim,T}(data, starttime, endtime, dt, par);
end

output = solvethis(setup_test_example2())


midnode = output.arrayoutputs["force"]
plot(midnode.u, midnode.fint, reuse=false, marker="o")=#