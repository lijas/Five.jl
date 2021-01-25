using Five

data = ProblemData(
    dim = 2,
    tend = 1.0
)

data.grid = generate_grid(Quadrilateral, (50,5), Vec((0.0,0.0)), Vec((1000.0, 50.0)))

addvertexset!(data.grid, "topright", (x) -> x[1] ≈ 1000.0 && x[2] ≈ 50.0)

material = MatHyperElasticPlastic(
    elastic_material = MatNeoHook(
        E = 200.0,
        ν = 0.3
    ),
    τ₀ = 1.0/10,
    H = 200.0/10
) |> PlaneStrainMaterial

material = MatNeoHook(
    E = 200.0,
    ν = 0.3
) |> PlaneStrainMaterial

#=material = MatLinearElastic(
    E = 1e5,
    nu = 0.3
)
=#
#=@addmat MatLinearElastic, "Steel", 1 begin
    E = 1e5,
    nu = 0.3   
end=#

con1 = Dirichlet(
    set = getfaceset(data.grid, "left"),
    func = (x,t) -> (0.0, 0.0),
    field = :u,
    dofs = [1,2]
)
push!(data.dirichlet, con1)

con1 = Dirichlet(
    set = getfaceset(data.grid, "right"),
    func = (x,t) -> (0.0,  t*10.0),
    field = :u,
    dofs = [1,2]
)
push!(data.dirichlet, con1)

part = Part{2,Float64}(
    element = Five.SolidElementQuad(),
    material = material,
    cellset = collect(1:getncells(data.grid))
)
push!(data.parts, part)

data.output[] = Output(
    interval = 0.1,
    runname = "Beamexample2",
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

force = PointForce(
    field = :u,
    comps = [2],
    set = getvertexset(data.grid, "topright"),
    func = (X,t) -> -1.0*t
)
#push!(data.external_forces, force)

solver = NewtonSolver(
    Δt0 = 0.1,
    Δt_max = 0.1,
    tol = 1e-5
)

state, data = build_problem(data)

solvethis(solver, state, data)
    
