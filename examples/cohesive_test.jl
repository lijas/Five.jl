using Five


function generate_test_grid(nelx, nely, L, h, a0)

    #
    grid1 = generate_grid(Ferrite.Quadrilateral,(nelx,nely),Vec((0.0,0.0)),Vec((L,h)))
    grid2 = generate_grid(Ferrite.Quadrilateral,(nelx,nely),Vec((0.0,h)),Vec((L,h*2)))
    grid = gridmerge(grid1,grid2)

    #
    addvertexset!(grid, "topvertices", (x)->  x[2] ≈ 2h)
    addfaceset!(grid, "topface", (x)->  x[2] ≈ 2h)
    addfaceset!(grid, "bot", (x)-> x[2] ≈ 0.0)

    #
    construct_interfacer_cells!(grid, "top1", "bottom2")

    #
    solid_cells = collect(1:nelx*nely*2)
    cz_cells = collect((1:nelx) .+ 2*nelx*nely)

    addcellset!(grid, "solid_cells", solid_cells)
    addcellset!(grid, "cz_cells", cz_cells)

    return grid
end

function construct_interfacer_cells!(grid, setname1::String, setname2::String)

    
    grid2_bottom_faceset = collect(getfaceset(grid, setname2))
    grid1_top_faceset = collect(getfaceset(grid, setname1))
    function sortby(f1)
        n1,n2 = Ferrite.faces(grid.cells[f1[1]])[f1[2]]
        return min(n1,n2)
    end
    function myless(f1,f2)
        n1,n2 = Ferrite.faces(grid.cells[f1[1]])[f1[2]]
        minA = min(grid.nodes[n1].x,grid.nodes[n1].x)

        n1,n2 = Ferrite.faces(grid.cells[f2[1]])[f2[2]]
        minB = min(grid.nodes[n1].x,grid.nodes[n1].x)
        return minA<minB
    end
    
    grid1_top_faceset = sort(grid1_top_faceset, lt = myless)
    grid2_bottom_faceset = sort(grid2_bottom_faceset, lt = myless)
    
    
    ncells = getncells(grid)
    for (i,faceind) in enumerate(grid1_top_faceset)
        botface_index = grid2_bottom_faceset[i]

        
        topface = Ferrite.faces(grid.cells[faceind[1]])[faceind[2]]
        botface = Ferrite.faces(grid.cells[botface_index[1]])[botface_index[2]]

        cell_nodes = [topface[2], topface[1], botface[1], botface[2]]
        new_cell = CohesiveCell{2,4,2}(Tuple(cell_nodes))
        push!(grid.cells, new_cell)
    end

end

#Dimension
const DIM = 2
const NELX = 2
const NELY = 2

const L = 2.0
const h = 2.0
const b = 2.0
const a0 = 46.9

data = ProblemData(
    dim = DIM,
    tend = .99,
    adaptive = true
)

#grid
data.grid = generate_test_grid(NELX, NELY, L, h, a0)

interfacematerial = 
MatCZBilinear(
    K    = 1.0e5,
    Gᴵ   = (0.5, 0.5, 0.5),
    τᴹᵃˣ = (50.0, 50.0, 50.0),
    η    = 1.6
) 

material = 
MatLinearElastic(
    E = 1e7,
    nu = 0.3
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
        nqp = 3,
        thickness = b
    ),
    material = interfacematerial,
    cellset = getcellset(data.grid, "cz_cells") |> collect
)
push!(data.parts, part)

#
dbc1 = Ferrite.Dirichlet(
    field = :u,
    set = getfaceset(data.grid, "bot"),
    func = (x,t)->[0.0, 0.0],
    dofs =  [1,2]
)
push!(data.dirichlet, dbc1)

dbc1 = Ferrite.Dirichlet(
    field = :u,
    set = getfaceset(data.grid, "topface"),
    func = (x,t)->[0.0, t*0.02],
    dofs =  [1,2]
)
push!(data.dirichlet, dbc1)

#
force = PointForce(
    field = :u,
    comps = [2],
    set = getvertexset(data.grid, "topvertices"),
    func = (X,t) -> 1.0/2 
)
#push!(data.external_forces, force)

#
data.output[] = Output(
    interval = 10.1,
    runname = "cohesiv_test",
    savepath = "."
)

#
output = OutputData(
    type = DofValueOutput(
        field = :u,
        dofs = 1:2
    ),
    interval = 0.0,
    set = getvertexset(data.grid, "topvertices")
)
data.outputdata["reactionforce"] = output

state, globaldata = build_problem(data)

solver = NewtonSolver(
    Δt0 = 0.001,
    Δt_max = 0.001,
)

#=solver = DissipationSolver(
    Δλ0          = 5.0,
    Δλ_max       = 10.0,
    Δλ_min       = 0.1,
    ΔL0          = 0.01,
    ΔL_min       = 1e-2,
    ΔL_max       = 0.2,
    sw2d         = 0.2,
    sw2i         = 1e-7,
    optitr       = 5,
    maxitr       = 10,
    maxsteps     = 200,
    λ_max        = 400.0,
    λ_min        = -100.0,
    tol          = 1e-8,
    max_residual = 1e5
)=#

output = solvethis(solver, state, globaldata)


d = [output.outputdata["reactionforce"].data[i].displacement for i in 1:length(output.outputdata["reactionforce"].data)]
f = [output.outputdata["reactionforce"].data[i].fint for i in 1:length(output.outputdata["reactionforce"].data)]
