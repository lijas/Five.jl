# # ENF example

using Five

include("enf_grid_generator.jl")

#Crack length
a0 = 30.9

#Celltype
CohesiveCellType = CohesiveCell{2,4,2} 
SolidCellType = Ferrite.Quadrilateral

#Dimension
DIM = 2
NELX = 176
NELY = 1

ORDERS = (2,2)

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
Five.TransverseIsotropicEngineeringConstants(
    E_L = 126e3, 
    E_T = 10e3, 
    G_LT = 8e3,
    ν_LT = 0.29, 
    ν_TT = 0.45
)

#
solidpart = Part(
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

#
part = Part(
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
    data.materialstates[cellid] = [Five.initial_material_state(material, Vec((1.0,0.0,0.0))) for i in 1:nqp]
end

#
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


#
force = PointForce(
    field = :u,
    comps = [2],
    set = [first(getvertexset(data.grid, "mid"))],
    func = (X,t) -> -1.0
)
push!(data.external_forces, force)

#
data.output[] = Output(
    interval = 0.0,
    runname = "enf_a0$(floor(Int,a0))_locdis_",
    savepath = "."
)

#
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
    type = Five.StressOutput(),
    func = mean,
)
Five.push_vtkoutput!(data.output[], vtkoutput)

vtkoutput = VTKCellOutput(
    type = MaterialStateOutput(
        field = :d,
        datatype = Float64
    ),
    func = mean,
)
Five.push_vtkoutput!(data.output[], vtkoutput)

vtkoutput = VTKNodeOutput(
    type = Five.StressOutput(),
)
Five.push_vtkoutput!(data.output[], vtkoutput)

state, globaldata = build_problem(data)

solver = LocalDissipationSolver(
    Δλ0          = 5.0,
    Δλ_max       = 40.0,
    Δλ_min       = 1e-7,
    ΔL0          = 2.5,
    ΔL_min       = 1e-2,
    ΔL_max       = 10.0,
    sw2d         = 0.5,
    sw2i         = 1e-7,
    optitr       = 8,
    maxitr       = 13,
    maxitr_first_step = 50,
    maxsteps     = 10,
    λ_max        = 1200.0,
    λ_min        = -100.0,
    tol          = 1e-8,
    max_residual = 1e5
)

output = solvethis(solver, state, globaldata)

d = [output.outputdata["reactionforce"].data[i].displacement[2] for i in 1:length(output.outputdata["reactionforce"].data)]
f = [output.outputdata["reactionforce"].data[i].fint[2] for i in 1:length(output.outputdata["reactionforce"].data)]