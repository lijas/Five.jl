# # Beam example

using Five

data = ProblemData(
    dim = 3,
    tend = 0.1
)

data.grid = generate_grid(Hexahedron, (50,20,20), Vec((0.0, 0.0, 0.0)), Vec((10.0, 4.0, 1.0)))

addvertexset!(data.grid, "topright", (x) -> x[1] == 10.0 && x[2] == 4.0  && x[3] == 1.0)

material = LinearElastic(
    E = 1e5,
    ν = 0.3
)

#=material = TransverselyIsotropic(;    
    EL = 100.0,   
    ET = 100.0,   

    ν_12 = 0.3, 
    G_12 = 100.0/(2(1.3)),
    ρ = 1.0,
    α =0.0
)=#

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
    func = (x,t) -> (0.0, 0.0, 0.0),
    field = :u,
    dofs = [1,2,3]
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
    element  = Five.LinearSolidElement{3,1,RefCube,Float64}(
        thickness = 1.0, 
        qr_order = 2,
        celltype = Hexahedron,
        #dimstate = PlaneStrain()
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
    comps = [3],
    set = getvertexset(data.grid, "topright"),
    func = (X,t) -> -10.0*t
)
push!(data.external_forces, force)

solver = NewtonSolver(
    Δt0 = 0.1,
    Δt_max = 0.1,
    tol = 1e-4,
    linearsolver = Five.LinearSolve.KrylovJL_CG(itmax = 100000, verbose=1, atol = 1e-8),
    #preconditioner = Five.Preconditioners.DiagonalPreconditioner
    preconditioner = Five.Preconditioners.AMGPreconditioner{Five.Preconditioners.RugeStuben}#Five.IncompleteLU.ilu#Five.IdentityProconditioner #Five.IncompleteLU.ilu
)

state, data = build_problem(data)

output = solvethis(solver, state, data)

d = output.outputdata["reactionforce"].data[end].displacement

using Test 
@show abs(d[1])                      #src
#@test norm(d) ≈ 0.44438080510979344 #src