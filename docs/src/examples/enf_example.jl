using Five
using Ferrite

include("enf_grid_generator.jl")

#Crack length
a0 = 16.9

#Celltype
CohesiveCellType = CohesiveCell{2,4,2}
SolidCellType = Ferrite.Quadrilateral

#Dimension
DIM = 2
NELX = 176
NELY = 1

L = 120.0
h = 2.0
b = 20.0

data = ProblemData(
    dim = DIM,
    tend = 1.0,
)

#grid
data.grid = generate_enf_grid(NELX, NELY, L, h, a0, SolidCellType, CohesiveCellType)

interfacematerial =
MatCZBilinear(
    K    = 1.0e5,
    Gᴵ   = (0.5, 0.5, 0.5),
    τᴹᵃˣ = (50.0, 50.0, 50.0),
    η    = 1.0
)

material =
MatTransvLinearElastic(
    E1 = 126.0e3,
    E2 = 10.0e3,
    ν_12 = 0.29,
    G_12 = 8.0e3,
    α = 0.0
)

solidpart = Part{2,Float64}(
    element  = Five.SolidElement{2,1,RefCube,Float64}(
        thickness = b,
        qr_order = 2,
        celltype = SolidCellType,
        dimstate = PlaneStrain()
    ),
    material = material,
    cellset  = getcellset(data.grid, "solid_cells") |> collect |> sort
)
push!(data.parts, solidpart)

part = Part{2,Float64}(
    element = CohesiveElement(
        order = 1,
        thickness = b,
        nqp = 2,
        celltype = CohesiveCellType
    ),
    material = interfacematerial,
    cellset = getcellset(data.grid, "cz_cells") |> collect |> sort
)
push!(data.parts, part)

#Change initial states of interface damage
for cellid in getcellset(data.grid, "precracked")
    nqp = 2 #Hardcoded
    data.materialstates[cellid] = [Five.initial_material_state(interfacematerial, 1.0) for i in 1:nqp]
end

#Add material directions for composite material
for cellid in solidpart.cellset
    nqp = 4 #hardcoded
    data.materialstates[cellid] = [Five.initial_material_state(material, Vec((1.0,0.0,0.0)) for i in 1:nqp)]
end

dbc1 = Ferrite.Dirichlet(
    field = :u,
    set = getvertexset(data.grid, "botleft"),
    func = (x,t)->[0.0, 0.0],
    dofs =  [1,2]
)
push!(data.dirichlet, dbc1)

dbc1 = Ferrite.Dirichlet(
    field = :u,
    set    = getvertexset(data.grid, "botright"),
    func   = (x,t)->[0.0],
    dofs  = [2]
)
push!(data.dirichlet, dbc1)

force = PointForce(
    field = :u,
    comps = [2],
    set = [first(getvertexset(data.grid, "mid"))],
    func = (X,t) -> -1.0
)
push!(data.external_forces, force)

data.output[] = Output(
    interval = 0.0,
    runname = "enf_a0$(floor(Int,a0))_locdis_",
    savepath = "."
)

output = OutputData(
    type = DofValueOutput(
        field = :u,
        dofs = 1:2
    ),
    interval = 0.0,
    set = Set([first(getvertexset(data.grid, "mid"))])
)
data.outputdata["reactionforce"] = output

vtkoutput = VTKNodeOutput(
    type = StressOutput(),
    func = mean,
)
Five.push_vtkoutput!(data.output[], vtkoutput)

vtkoutput = VTKCellOutput(
    type = Five.StressOutput(),
)
Five.push_vtkoutput!(data.output[], vtkoutput)

state, globaldata = build_problem(data)

solver = LocalDissipationSolver(
    Δλ0          = 5.0,
    Δλ_max       = 10.0,
    Δλ_min       = 1e-7,
    ΔL0          = 2.5,
    ΔL_min       = 1e-2,
    ΔL_max       = 5.0,
    sw2d         = 0.2,
    sw2i         = 1e-7,
    optitr       = 8,
    maxitr       = 13,
    maxitr_first_step = 50,
    maxsteps     = 200,
    λ_max        = 1200.0,
    λ_min        = -100.0,
    tol          = 1e-5,
    max_residual = 1e5
)

output = solvethis(solver, state, globaldata)

d = [output.outputdata["reactionforce"].data[i].displacement for i in 1:length(output.outputdata["reactionforce"].data)]
f = [output.outputdata["reactionforce"].data[i].fint for i in 1:length(output.outputdata["reactionforce"].data)]

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

