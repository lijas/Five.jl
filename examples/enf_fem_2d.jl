using Five


function generate_enf_grid(nelx, nely, L, h, a0)

    #
    grid1 = generate_grid(JuAFEM.Quadrilateral,(nelx,nely),Vec((0.0,0.0)),Vec((L,h)))
    grid2 = generate_grid(JuAFEM.Quadrilateral,(nelx,nely),Vec((0.0,h)),Vec((L,h*2)))
    grid = gridmerge(grid1,grid2)

    #
    addvertexset!(grid, "mid", (x)-> x[1] ≈ L/2 && x[2] ≈ h*2)
    @assert(length(getvertexset(grid, "mid")) == 2)
    addvertexset!(grid, "botleft", (x)-> x[1] ≈ 0.0 && x[2] ≈ 0.0)
    addvertexset!(grid, "botright", (x)-> x[1] ≈ L && x[2] ≈ 0.0)

    #
    construct_interfacer_cells!(grid, "top1", "bottom2")

    #
    solid_cells = collect(1:nelx*nely*2)
    cz_cells = collect((1:nelx) .+ 2*nelx*nely)

    addcellset!(grid, "solid_cells", solid_cells)
    addcellset!(grid, "cz_cells", cz_cells)

    addcellset!(grid, "temp_precraced", (x)-> x[1]>L-a0)
    precracked_cells = setdiff(getcellset(grid, "temp_precraced"), solid_cells)

    addcellset!(grid, "precracked", precracked_cells)

    return grid
end

function construct_interfacer_cells!(grid, setname1::String, setname2::String)

    
    grid2_bottom_faceset = collect(getfaceset(grid, setname2))
    grid1_top_faceset = collect(getfaceset(grid, setname1))
    function sortby(f1)
        n1,n2 = JuAFEM.faces(grid.cells[f1[1]])[f1[2]]
        return min(n1,n2)
    end
    function myless(f1,f2)
        n1,n2 = JuAFEM.faces(grid.cells[f1[1]])[f1[2]]
        minA = min(grid.nodes[n1].x,grid.nodes[n1].x)

        n1,n2 = JuAFEM.faces(grid.cells[f2[1]])[f2[2]]
        minB = min(grid.nodes[n1].x,grid.nodes[n1].x)
        return minA<minB
    end
    
    grid1_top_faceset = sort(grid1_top_faceset, lt = myless)
    grid2_bottom_faceset = sort(grid2_bottom_faceset, lt = myless)
    
    
    ncells = getncells(grid)
    for (i,faceind) in enumerate(grid1_top_faceset)
        botface_index = grid2_bottom_faceset[i]

        
        topface = JuAFEM.faces(grid.cells[faceind[1]])[faceind[2]]
        botface = JuAFEM.faces(grid.cells[botface_index[1]])[botface_index[2]]

        cell_nodes = [topface[2], topface[1], botface[1], botface[2]]
        new_cell = Cell{2,4,1}(Tuple(cell_nodes))
        push!(grid.cells, new_cell)
    end

end

#Dimension
 DIM = 2
 NELX = 176
 NELY = 1

 ORDERS = (2,2)

 L = 120.0
 h = 2.0
 b = 20.0
 a0 = 46.9

angles = deg2rad.([0.0, 0.0])
nlayers = length(angles)
ninterfaces = nlayers - 1

data = ProblemData(
    dim = DIM,
    tend = 1.0,
)

#grid
data.grid = generate_enf_grid(NELX, NELY, L, h, a0)

#interfacematerial = IGAShell.MatCohesive{dim}(λ_0,λ_f,τ,K)
interfacematerial = 
MatCZBilinear(
    K    = 1.0e5,
    Gᴵ   = (0.5, 0.5, 0.5),
    τᴹᵃˣ = (50.0, 50.0, 50.0),
    η    = 1.0
) 

#interfacematerial = MatMixedCohesive(1.0e5, (0.5, 0.5), (50.0,50.0), 1.0) 

#=
interfacematerial = 
MatCZBilinear2(
    K    = 1.0e5,
    Gᴵ   = (0.5, 0.5, 0.5),
    τᴹᵃˣ = (50.0, 50.0, 50.0),
    η    = 1.0
) =#

material = 
MatTransvLinearElastic(
    E1 = 126.0e3,
    E2 = 10.0e3,
    ν_12 = 0.29, 
    G_12 = 8.0e3, 
    α = 0.0
) |> PlaneStrainMaterial

#
part = Part{2,Float64}(
    element  = SolidElementQuad(thickness = b),
    material = material,
    cellset  = getcellset(data.grid, "solid_cells") |> collect
)
push!(data.parts, part)

#
part = Part{2,Float64}(
    element = CohesiveElement{DIM}(
        order = 1,
        thickness = b
    ),
    material = interfacematerial,
    cellset = getcellset(data.grid, "cz_cells") |> collect
)
push!(data.parts, part)

#Change initial states
for cellid in getcellset(data.grid, "precracked")
    data.materialstates[cellid] = [Five.getmaterialstate(interfacematerial, 1.0) for i in 1:2]
end

#
dbc1 = JuAFEM.Dirichlet(
    field = :u,
    set = getvertexset(data.grid, "botleft"),
    func = (x,t)->[0.0, 0.0],
    dofs =  [1,2]
)
push!(data.dirichlet, dbc1)

dbc1 = JuAFEM.Dirichlet(
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
    runname = "enf_example",
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

state, globaldata = build_problem(data)

solver = DissipationSolver(
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
    λ_max        = 400.0,
    λ_min        = -100.0,
    tol          = 1e-3,
    max_residual = 1e5
)

output = solvethis(solver, state, globaldata)

d = [output.outputdata["reactionforce"].data[i].displacement for i in 1:length(output.outputdata["reactionforce"].data)]
f = [output.outputdata["reactionforce"].data[i].fint for i in 1:length(output.outputdata["reactionforce"].data)]
