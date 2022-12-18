using Five


#NELX = 51
#NELY = 51

b = 0.001
E = 6e11
nu = 0.22
#elsize = L/NELY

data = ProblemData(
    dim = 2,
    tend = 0.010,
)

#grid
include("phase_grid.jl")

h1 = Vec((0.02, 0.02))
h2 = Vec((0.02, 0.1))

data.grid = read_grid(joinpath(@__DIR__, "paperb_phase_mesh.mphtxt"))#generate_grid(Quadrilateral, (NELX,NELY), Vec((0.0,0.0)),Vec((L,L)))
#addvertexset!(data.grid, "leftbottom", x -> x[2] ≈ 0.0 && x[1] ≈ 0.0)
#addvertexset!(data.grid, "top", x -> x[2] ≈ 1.0)
addfaceset!(data.grid, "bottomhole", x ->  norm(x-h1) <= .010)
addfaceset!(data.grid, "tophole", x ->  norm(x-h2) <= .010)
addvertexset!(data.grid, "tophole", x ->  norm(x-h2) <= .010)

material = Five.PhaseFieldSpectralSplit(
    E = E,
    ν = nu,    
    Gc = 2280.0,
    lc = 0.25e-3*2,
)

#
element = PhaseFieldElement{2,1,RefCube,Float64}(
    thickness = 0.001,     
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
    set = getfaceset(data.grid, "bottomhole"),
    func = (x,t)->[0.0, 0.0],
    dofs =  [1, 2]
)
push!(data.dirichlet, dbc1)

dbc1 = Ferrite.Dirichlet(
    field = :u,
    set = getfaceset(data.grid, "tophole"),
    func = (x,t)->[t*1.0],
    dofs =  [2]
)
#push!(data.dirichlet, dbc1)

#
force = PointForce(
    field = :u,
    comps = [2],
    set = getvertexset(data.grid, "tophole"),
    func = (X,t) -> [1.0]
)
push!(data.external_forces, force)


#
data.output[] = Output(
    interval = -1.0,
    runname = "phase_comsol",
    savepath = "."
)

#
output = OutputData(
    type = DofValueOutput(
        field = :u,
        dofs = [2]
    ),
    interval = -1.0,
    set = getfaceset(data.grid, "tophole")
)
data.outputdata["reactionforce"] = output

state, globaldata = build_problem(data)

solver = NewtonSolver(
    Δt0 = 0.0000001,
    Δt_max = 0.0001,
    tol    = 1e-4,
)


solver = LocalDissipationSolver(
    Δλ0          = 1e-7,
    Δλ_max       = 10.0,
    Δλ_min       = 1e-10,
    #Δλ0          = 1.0,
    #Δλ_max       = 10.0,
    #Δλ_min       = 1e-11,
    ΔL0          = 2.5,
    ΔL_min       = 1e-9,
    ΔL_max       = 10.0,
    sw2d         = 1e-6,
    sw2i         = -Inf,
    optitr       = 7,
    maxitr       = 12,
    maxitr_first_step = 50,
    maxsteps     = 100,
    λ_max        = 100.108,
    #λ_max        = 100.0,
    λ_min        = -2000.0,
    tol          = 1e-8,
    max_residual = 1e8,
    #finish_criterion = Five.finish_criterion1
)


output = solvethis(solver, state, globaldata)

d = [output.outputdata["reactionforce"].data[i].displacement for i in 1:length(output.outputdata["reactionforce"].data)]
f = [output.outputdata["reactionforce"].data[i].fint for i in 1:length(output.outputdata["reactionforce"].data)]

#ig = plot(xlabel = "Displacement", ylabel = "Force", legend = false)
#plot!(fig, d, f, mark = :o)
#plot!(fig, d1, f1, label = "Elastic energy", mark = :o)
