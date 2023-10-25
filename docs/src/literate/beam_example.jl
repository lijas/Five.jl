# # Beam example

using Five

data = ProblemData(
    runname = "Beamexample",
    dim = 2,
    tend = 1.0,
    vtk_output_interval = 0.1,
    vtkoutputtype = Five.FerriteVTKOutput(),
)

data.grid = generate_grid(Quadrilateral, (10,5), Vec((0.0, 0.0)), Vec((10.0, 1.0)))
addvertexset!(data.grid, "topright", (x) -> x[1] == 10.0 && x[2] == 1.0)

vtkoutput = VTKNodeOutput(
    type = StressOutput()
)
push!(data.vtk_output, vtkoutput)

vtkoutput = VTKNodeOutput(
    type = MaterialStateOutput(
        field    = :κ,
        datatype = Float64
    )
)
push!(data.vtk_output, vtkoutput)

output = OutputData(
    type = DofValueOutput(
        field = :u,
        dofs = [2]
    ),
    interval = 0.1,
    set = getvertexset(data.grid, "topright")
)
data.outputdata["reactionforce"] = output

material = LinearElastic(;    
    E = 100.0,   
    ν = 0.3,   
)

material = Plastic(E=200e3, ν=0.3, σ_y=200., H=50., r=0.5, κ_∞=13., α_∞=13.)

con1 = Dirichlet(
    set = getfaceset(data.grid, "left"),
    func = (x,t) -> (0.0, 0.0),
    field = :u,
    dofs = [1,2]
)
push!(data.constraints_ferrite, con1)

part = Part(
    element  = Five.LinearSolidElement(
        celltype = Quadrilateral,
        qr_order = 2,
        thickness = 1.0, 
        dimstate = PlaneStrain()
    ),
    material = material,
    cellset = collect(1:getncells(data.grid))
)
push!(data.parts, part)

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
    #linearsolver = Five.LinearSolve.KrylovJL_CG(itmax = 1000, verbose=0, atol = 1e-15),
    #preconditioner = Five.Preconditioners.AMGPreconditioner{Five.Preconditioners.RugeStuben}#Five.IncompleteLU.ilu#Five.IdentityProconditioner #Five.IncompleteLU.ilu
)

state, data = build_problem(data)

output = solvethis(solver, state, data)

d = output.outputdata["reactionforce"].data[end].displacement

using Test 
@show abs(d[1])                      #src
#@test norm(d) ≈ 0.44438080510979344 #src