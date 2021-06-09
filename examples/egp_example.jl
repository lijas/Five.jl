using Five

tmax = 1000.0
tsim = 1000.0
data = ProblemData(
    dim = 3,
    tend = tsim
)

data.grid = generate_grid(Hexahedron, (1,1,1), Vec((0.0, 0.0, 0.0)), Vec((1.0, 1.0, 1.0)))

addvertexset!(data.grid, "x100", (x) -> x[1] == 1.0 && x[2] == 0.0 && x[3] == 0.0)
addvertexset!(data.grid, "x000", (x) -> x[1] == 0.0 && x[2] == 0.0 && x[3] == 0.0)
addvertexset!(data.grid, "x111", (x) -> x[1] == 1.0 && x[2] == 1.0 && x[3] == 1.0)

material = 
MatEGP(
    E = 1.0e6,
    ν = 0.25,
    
    ACTION = 1,
    NOP = 1,
    NAM = 1,
    NBM = 0,
    NGM = 0,
    NGENS = 6,
    PRESMET = 1,
    HARDMET = 1,
    VISCHARD = 0,
    STIFFNESS = 0,
    VISCDEF = 0,

    STANDPROP = [NaN, 273, 273, 26.0, 0.0E0, 0.0E0, 3750, 1.0],

    PROPS_PROC1_STD = [2.89e5, 0.0, 0.0 , 0.0 , 26.5, 0.964, 50.0, -5.0, 5.38452923e-27, 0.08],
    PROPS_PROC1_G = [321.0],
    PROPS_PROC1_GH0R = [2.1e11]
) 

#=material = MatLinearElastic(
    E = 1e4,
    nu = 0.35
) =#

con1 = Dirichlet(
    set = getfaceset(data.grid, "front"),
    func = (x,t) -> (0.0),
    field = :u,
    dofs = [2]
)
push!(data.dirichlet, con1)

con1 = Dirichlet(
    set = getvertexset(data.grid, "x000"),
    func = (x,t) -> (0.0,0.0),
    field = :u,
    dofs = [1,3]
)
push!(data.dirichlet, con1)

con1 = Dirichlet(
    set = getvertexset(data.grid, "x100"),
    func = (x,t) -> (0.0),
    field = :u,
    dofs = [3]
)
push!(data.dirichlet, con1)

loadf(t) = begin
    if t < tmax/2
        return (-t/tmax) * -0.8
    else
        return -0.5 * -0.8#(t-tmax)/tmax * -0.8
    end
end

con1 = Dirichlet(
    set = getfaceset(data.grid, "back"),
    func = (x,t) -> exp(-0.001*t)-1.0,
    field = :u,
    dofs = [2]
)
push!(data.dirichlet, con1)

part = Part{3,Float64}(
    element = SolidElement{3,1,RefCube,Float64}(
        celltype = Ferrite.Hexahedron,
        qr_order = 2,
        total_lagrangian = false
    ),
    material = material,
    cellset = collect(1:getncells(data.grid))
)
push!(data.parts, part)

data.output[] = Output(
    interval = tmax/10,
    runname = "zegptest",
    savepath = "."
)

vtkoutput = VTKCellOutput(
    type = MaterialStateOutput(
        field = :εᵖ
    ),
    func = (x)-> mean(x),
)
Five.push_vtkoutput!(data.output[], vtkoutput)

vtkoutput = VTKCellOutput(
    type = MaterialStateOutput(
        field = :σ
    ),
    func = (x)-> mean(x),
)
Five.push_vtkoutput!(data.output[], vtkoutput)

output = OutputData(
    type = DofValueOutput(
        field = :u,
        dofs = [2]
    ),
    interval = 1.0,
    set = getvertexset(data.grid, "x111")
)
data.outputdata["reactionforce"] = output

solver = NewtonSolver(
    Δt0 = 10.0,
    Δt_min = 0.01,
    Δt_max = 10.0,
    tol = 1e-4,
    maxitr = 20,
    optitr = 15
)

state, data = build_problem(data)

output = solvethis(solver, state, data)
  
d = [output.outputdata["reactionforce"].data[i].displacement for i in 1:length(output.outputdata["reactionforce"].data)]
f = [output.outputdata["reactionforce"].data[i].fint for i in 1:length(output.outputdata["reactionforce"].data)]

using Plots; plotly()
fig = plot(abs.(d),abs.(f), mark=:o, xlim = [0.0, 1.0], ylim = [0.0,100.0], label="julia UpdateLag")
plot!(fig, abs.(dtl),abs.(ftl) , mark=:cross, xlim = [0.0, 1.0], ylim = [0.0,100.0], label="julia TotalLag")
plot!(fig, abs.(res[:,1]), abs.(res[:,2]) , mark=:square, label = "pyfem")
