using Five

data = ProblemData(
    dim = 3,
    tend = 1.0
)

data.grid = generate_grid(Hexahedron, (10,5,5), Vec((0.0, 0.0, 0.0)), Vec((10.0, 1.0, 1.0)))

addvertexset!(data.grid, "topright", (x) -> x[1] == 10.0 && x[2] == 1.0 && x[3] == 1.0)

E = 1e5
ν = 0.3

material = MaterialModels.NeoHook(
    λ = (E*ν)/((1+ν)*(1-2ν)),
    μ = E/(2*(1+ν))
)


material = LinearElastic(
    E = E,
    ν = ν
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
    func = (x,t) -> (0.0, 0.0,0.0),
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

part = Part{3,Float64}(
    element = SolidElement{3,1,RefCube,Float64}(
        celltype = Ferrite.Hexahedron,
        qr_order = 2,
        dimension = Five.ThreeD()
    ),
    material = material,
    cellset = collect(1:getncells(data.grid)),
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
        dofs = [3]
    ),
    interval = 0.1,
    set = getvertexset(data.grid, "topright")
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
    comps = [3],
    set = getvertexset(data.grid, "topright"),
    func = (X,t) -> -1.0*t
)
push!(data.external_forces, force)

solver = NewtonSolver(
    Δt0 = 0.1,
    Δt_max = 0.1,
    tol = 1e-5
)

state, data = build_problem(data)

output = solvethis(solver, state, data)

d = output.outputdata["reactionforce"].data[end].displacement

using Test                         #src
@test abs(d) ≈ 0.44438080510979344 #src