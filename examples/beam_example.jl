using Five

data = ProblemData(
    dim = 2,
    tend = 7.0
)

data.grid = generate_grid(Quadrilateral, (10,1), Vec((0.0, 0.0)), Vec((100.0, 10.0)))

addvertexset!(data.grid, "topright", (x) -> x[1] == 100.0 && x[2] == 10.0)

material = MatLinearElastic(
    density = 7e-6,
    E = 100.0,
    nu = 0.3
) |> PlaneStrainMaterial 

#=material = MatHyperElasticPlastic(
    elastic_material = MatNeoHook(
        E = 1.0e5,
        ν = 0.3
    ),
    τ₀ = 400.0,
    H = 1.0e5/20
) |> PlaneStrainMaterial=#

con1 = Dirichlet(
    set = getfaceset(data.grid, "right"),
    func = (x,t) -> (0.0, 0.0),
    field = :u,
    dofs = [1,2]
)
push!(data.dirichlet, con1)

con1 = Dirichlet(
    set = getfaceset(data.grid, "left"),
    func = (x,t) -> (0.0, 0.1*t^2),
    field = :u,
    dofs = [1,2]
)
push!(data.dirichlet, con1)

part = Part{2,Float64}(
    element = SolidElement{2,2,RefCube,Float64}(
        celltype = Ferrite.Quadrilateral,
        qr_order = 4
    ),
    material = material,
    cellset = collect(1:getncells(data.grid))
)
push!(data.parts, part)

data.output[] = Output(
    interval = 0.1,
    runname = "Beamexample",
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

#=
vtkoutput = VTKCellOutput(
    type = MaterialStateOutput(
        field = :ϵᵖ
    ),
    func = mean,
)
Five.push_vtkoutput!(data.output[], vtkoutput)

vtkoutput = VTKNodeOutput(
    type = MaterialStateOutput(
        field = :ϵᵖ
    ),
    func = mean,
)
Five.push_vtkoutput!(data.output[], vtkoutput)
=#

#=force = PointForce(
    field = :u,
    comps = [2],
    set = getvertexset(data.grid, "topright"),
    func = (X,t) -> -10.0*t
)
push!(data.external_forces, force)=#

solver = ExplicitSolver(
    Δt0 = 0.0001,
    #Δt_max = 0.1,
)

state, data = build_problem(data)

solvethis(solver, state, data)
    
