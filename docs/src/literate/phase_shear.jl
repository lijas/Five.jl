using Five

L = 1.0
a0 = L/2
E = 210.0e3
nu = 0.3

data = ProblemData(
    dim = 2,
    tend = 0.010,
)

#grid
include("phase_grid.jl")

data.grid = read_grid(joinpath(@__DIR__, "shear.mphtxt"))
data.grid.nodes .= [Node(Vec{2,Float64}(((node.x[1] + 0.0005)*1000, (node.x[2] + 0.0005)*1000))) for node in data.grid.nodes]
addvertexset!(data.grid, "topright", x -> x[2] ≈ L && x[1] ≈ L)
addvertexset!(data.grid, "top", x -> x[2] ≈ L)
addfaceset!(data.grid, "top", x -> x[2] ≈ L)
addfaceset!(data.grid, "bottom", x ->  x[2] ≈ 0.0)
addfaceset!(data.grid, "right", x -> x[1] ≈ L)
addfaceset!(data.grid, "left", x -> x[1] ≈ 0.0)

material = Five.PhaseFieldSpectralSplit(
    E = E,
    ν = nu,    
    Gc = 2.7,
    lc = 0.015, # ellength = 0.007323999999999997
)

#
element = PhaseFieldElement{2,1,RefCube,Float64}(
    thickness = 1.0,     
    qr_order = 2, 
    celltype = Quadrilateral,
    dimstate = MaterialModels.PlaneStrain()
)

#=element = Five.LinearSolidElement{2,1,RefCube,Float64}(
    thickness = 1.0,
    celltype = Quadrilateral,
    qr_order = 2,
    dimstate = PlaneStrain()
)=#

#
part = Part{Float64}(
    element = element,
    material = material,
    cellset  = 1:getncells(data.grid)
)
push!(data.parts, part)


#Change initial states
#for cellid in getcellset(data.grid, "precracked")
#    data.elementstates[cellid] = Five.PhaseFieldElementState(element, 1000.0)
#end

#
dbc1 = Ferrite.Dirichlet(
    field = :u,
    set = getfaceset(data.grid, "bottom"),
    func = (x,t)->[0.0, 0.0],
    dofs =  [1, 2]
)
push!(data.dirichlet, dbc1)


fc = Five.FollowerConstraint(
    field     = :u,
    component = 2,
    faces     = getfaceset(data.grid, "top"),
    mastervertex = first(getvertexset(data.grid, "topright")),
)
push!(data.dirichlet, fc)


dbc1 = Ferrite.Dirichlet(
    field = :u,
    set = getfaceset(data.grid, "left"),
    func = (x,t)->[0.0],
    dofs =  [2]
)
push!(data.dirichlet, dbc1)

dbc1 = Ferrite.Dirichlet(
    field = :u,
    set = getfaceset(data.grid, "right"),
    func = (x,t)->[0.0],
    dofs =  [2]
)
push!(data.dirichlet, dbc1)

force = PointForce(
    field = :u,
    comps = [1],
    set = getvertexset(data.grid, "top"),
    func = (X,t) -> [1.0]
)
push!(data.external_forces, force)

#
data.output[] = Output(
    interval = -1.0,
    runname = "phase_shear",
    savepath = "."
)

#
output = OutputData(
    type = DofValueOutput(
        field = :u,
        dofs = [1]
    ),
    interval = -1.0,
    set = getfaceset(data.grid, "top")
)
data.outputdata["reactionforce"] = output

state, globaldata = build_problem(data)


solver = LocalDissipationSolver(
    Δλ0          = 1e-5,
    Δλ_max       = 1e10,
    Δλ_min       = 1e-10,
    #Δλ0          = 1.0,
    #Δλ_max       = 10.0,
    #Δλ_min       = 1e-11,
    ΔL0          = 2.5,
    ΔL_min       = 1e-9,
    ΔL_max       = 1e10,
    sw2d         = 0.001,
    sw2i         = -Inf,
    optitr       = 7,
    maxitr       = 12,
    maxitr_first_step = 50,
    maxsteps     = 2,
    λ_max        = 1e10,
    #λ_max        = 100.0,
    λ_min        = -1e-10,
    tol          = 1e-6,
    max_residual = 1e10,
    #finish_criterion = Five.finish_criterion1
)


output = solvethis(solver, state, globaldata)

#d = [output.outputdata["reactionforce"].data[i].displacement for i in 1:length(output.outputdata["reactionforce"].data)]
#f = [output.outputdata["reactionforce"].data[i].fint for i in 1:length(output.outputdata["reactionforce"].data)]

#fig = plot(xlabel = "Displacement", ylabel = "Force", legend = false)
#plot!(fig, d, f, mark = :o)
#plot!(fig, d1, f1, label = "Elastic energy", mark = :o)