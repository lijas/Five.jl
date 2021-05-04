using Five


function generate_2dbuckle_grid(nelx, nely, L, h, d1)

    #
    h2 = h
    grid1 = generate_grid(Ferrite.Quadrilateral,(nelx, nely), Vec((0.0,0.0)), Vec((L,h2)))
    grid2 = generate_grid(Ferrite.Quadrilateral,(nelx, nely), Vec((0.0,d1+h2)), Vec((L,h+d1+h2)))
    grid3 = generate_grid(Ferrite.Quadrilateral,(nely÷3, round(Int, nelx/L *d1)), Vec((L - (L/nelx)*(nely÷3), h2)), Vec((L, h2+d1)))
    grid  = gridmerge(grid1, grid2, grid3)

    #Filter and merge coincident nodes
    nodepairs = Dict{Int,Int}()
    for i in eachindex(grid.nodes)
        for j in i:length(grid.nodes)
            i==j && continue

            if grid.nodes[i].x == grid.nodes[j].x
                nodepairs[i] = j
            end
        end
    end

    CellType = Ferrite.Quadrilateral
    for (cellid, cell) in enumerate(grid.cells)
        nodes_to_switch = []
        for (i, nodeid) in enumerate(cell.nodes)
            if nodeid in keys(nodepairs)
                newnodeid = nodepairs[nodeid]
                push!(nodes_to_switch, i=>newnodeid)
            end
        end

        #Found nodes taht should be switched
        if length(nodes_to_switch) != 0
            cellnodes = collect(cell.nodes)
            for pair in nodes_to_switch
                cellnodes[pair[1]] = pair[2]
            end
            newcell = CellType(Tuple(cellnodes))
            grid.cells[cellid] = newcell
        end

    end

    #disturb some nodes
    x0 = Vec( (L/2, h2/2) )
    r = h*1
    #addnodeset!(grid, "move", (x) -> norm(x0 - x) < r )
    addnodeset!(grid, "move", (x) -> x0[1]-r < x[1] < x0[1]+r && x[2]< h*1.1)
    for nodeid in getnodeset(grid, "move")
        x, y = grid.nodes[nodeid].x |> collect
        
        dir = x0 - grid.nodes[nodeid].x
        len = norm(dir)
        if len == 0
            continue
        end
        scl = (x0[2]-y)/(h2/2)
        newpos = Vec((x,y + scl*h/5))#grid.nodes[nodeid].x + (dir/len)*(len*h/(1))

        grid.nodes[nodeid] = Node(newpos)
    end

    #
    addfaceset!(grid, "left", (x)-> x[1] ≈ 0.0)
    addvertexset!(grid, "mid", (x)-> x[1] ≈ L && x[2] ≈ 0.0)

    #
    solid_cells = 1:getncells(grid) |> collect

    addcellset!(grid, "solid_cells", solid_cells)

    return grid
end


#Dimension
 DIM = 2
 NELX = 100
 NELY = 6

 L = 10.0
 h = 0.4/2
 b = 1.0
 d1 = 0.5

data = ProblemData(
    dim = DIM,
    tend = 1.0,
)

#grid
data.grid = generate_2dbuckle_grid(NELX, NELY, L, h, d1)

material = 
MatLinearElastic(
    E = 126.0e3,
    nu = 0.0,
    plane_stress = false
) 

material = MatHyperElasticPlastic(
    elastic_material = MatNeoHook(
        E = 1.0e5,
        ν = 0.3
    ),
    τ₀ = 400.0,
    H = 1.0e5/20
) |> PlaneStrainMaterial

#
part = Part{2,Float64}(
    element  = SolidElementQuad(thickness = b, qr_order = 2, celltype = Quadrilateral),
    material = material,
    cellset  = getcellset(data.grid, "solid_cells") |> collect
)
push!(data.parts, part)

#
dbc1 = Ferrite.Dirichlet(
    field = :u,
    set = getfaceset(data.grid, "left"),
    func = (x,t)->[0.0, 0.0],
    dofs =  [1,2]
)
push!(data.dirichlet, dbc1)

#
force = PointForce(
    field = :u,
    comps = [1],
    set = [first(getvertexset(data.grid, "mid"))],
    func = (X,t) -> -1.0
)
push!(data.external_forces, force)

#
data.output[] = Output(
    interval = 0.0,
    runname = "plastbuckle",
    savepath = "."
)

#
output = OutputData(
    type = DofValueOutput(
        field = :u,
        dofs = [1],
    ),
    interval = 0.0,
    set = Set([first(getvertexset(data.grid, "mid"))])
)
data.outputdata["reactionforce"] = output

vtkoutput = VTKNodeOutput(
    type = MaterialStateOutput(
        field = :ϵᵖ
    ),
    func = mean,
)
Five.push_vtkoutput!(data.output[], vtkoutput)

vtkoutput = VTKCellOutput(
    type = MaterialStateOutput(
        field = :ϵᵖ
    ),
    func = mean,
)
Five.push_vtkoutput!(data.output[], vtkoutput)

state, globaldata = build_problem(data)

solver = LocalDissipationSolver(
    Δλ0          = 0.1,
    Δλ_max       = 10.0,
    Δλ_min       = 1e-7,
    ΔL0          = 2.5,
    ΔL_min       = 1e-4,
    ΔL_max       = 7.0,
    sw2d         = 1e-5,
    sw2i         = 1e-7,
    optitr       = 5,
    maxitr       = 10,
    maxitr_first_step = 50,
    maxsteps     = 200,
    λ_max        = 20.0,
    λ_min        = 5.0,
    tol          = 1e-4,
    max_residual = 1e5,
    finish_criterion = Five.finish_criterion1
)

#=solver = ArcLengthSolver(
    Δλ0          = 0.1,

    λ_max        = 20.0,
    λ_min        = 5.0,

    ΔL_min       = 1e-4,
    ΔL_max       = 7.0,

    ψ = 0.0,
    maxsteps = 100,

    #Newton
    tol = 1.0e-5,
    max_residual = 1e5,
    optitr = 5,
    maxitr = 10,
    maxitr_first_step = 50,
    finish_criterion = Five.finish_criterion1
)=#

output = solvethis(solver, state, globaldata)

d = [output.outputdata["reactionforce"].data[i].displacement for i in 1:length(output.outputdata["reactionforce"].data)]
f = [output.outputdata["reactionforce"].data[i].fint for i in 1:length(output.outputdata["reactionforce"].data)]
