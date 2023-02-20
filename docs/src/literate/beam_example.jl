# # Beam example

using Five

data = ProblemData(
    dim = 2,
    tend = 1.0
)

data.grid = generate_grid(Quadrilateral, (10,5), Vec((0.0, 0.0)), Vec((10.0, 1.0)))

addvertexset!(data.grid, "topright", (x) -> x[1] == 10.0 && x[2] == 1.0)

#=material = MatLinearElastic(
    E = 1e5,
    nu = 0.3
)=#

material = MatTransvLinearElastic(;    
    E1 = 100.0,   
    E2 = 100.0,   

    ν_12 = 0.3, 
    G_12 = 100.0/(2(1.3)),
    ρ = 1.0,
    α =0.0
)

#= material = MatHyperElasticPlastic(
    elastic_material = MatNeoHook(
        E = 1.0e5,
        ν = 0.3
    ),
    τ₀ = 400.0,
    H = 1.0e5/20
) |> PlaneStrainMaterial =#

con1 = Dirichlet(
    set = getfaceset(data.grid, "left"),
    func = (x,t) -> (0.0, 0.0),
    field = :u,
    dofs = [1,2]
)
push!(data.dirichlet, con1)

#=
con1 = Dirichlet(
    set = getfaceset(data.grid, "right"),
    func = (x,t) -> (t*1.0),
    field = :u,
    dofs = [2]
)
push!(data.dirichlet, con1)=#

part = Part(
    element  = Five.LinearSolidElement{2,1,RefCube,Float64}(
        thickness = 1.0, 
        qr_order = 2,
        celltype = Quadrilateral,
        dimstate = PlaneStrain()
    ),
    material = material,
    cellset = collect(1:getncells(data.grid))
)
push!(data.parts, part)

data.output[] = Output(
    interval = -0.1,
    runname = "Beamexample",
    savepath = "."
)

output = OutputData(
    type = DofValueOutput(
        field = :u,
        dofs = [2]
    ),
    interval = 0.1,
    set = getvertexset(data.grid, "topright")
)
data.outputdata["reactionforce"] = output

vtkoutput = VTKNodeOutput(
    type = MaterialStateOutput(
        field = :σ,
        datatype = SymmetricTensor{2,3,Float64,6}
    ),
    func = mean,
)
Five.push_vtkoutput!(data.output[], vtkoutput)

#=vtkoutput = VTKNodeOutput(
    type = MaterialStateOutput(
        field = :ϵᵖ
    ),
    func = mean,
)
Five.push_vtkoutput!(data.output[], vtkoutput)=#


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
    linearsolver = Five.LinearSolve.KrylovJL_CG(itmax = 1000, verbose=0, atol = 1e-15),
    preconditioner = Five.Preconditioners.AMGPreconditioner{Five.Preconditioners.RugeStuben}#Five.IncompleteLU.ilu#Five.IdentityProconditioner #Five.IncompleteLU.ilu
)

state, data = build_problem(data)

output = solvethis(solver, state, data)

d = output.outputdata["reactionforce"].data[end].displacement

using Test 
@show abs(d[1])                      #src
#@test norm(d) ≈ 0.44438080510979344 #src