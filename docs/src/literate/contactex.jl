# # Beam example

using Five


function generate_my_grid2()
    
    grid1 = generate_grid(Quadrilateral, (5,2), Vec((-20.0, 1.0+0.00000001)), Vec((20.0, 2.0))) 
    #addfaceset!(grid2, "bottom", (x)->x[2]≈ 4.0)
    #addfaceset!(grid2, "top", (x)->x[2]≈5.0)
    #addfaceset!(grid2, "left", (x)->x[1]≈-20.0)
    #addfaceset!(grid2, "right", (x)->x[1]≈20.0)
    addnodeset!(grid1, "bottom", (x)->x[2]≈1.0+0.00000001)
    addcellset!(grid1, "bottomcells", (x)->true)

    grid2 = generate_grid(Quadrilateral, (2,2), Vec((-5.0,0.0)), Vec((5.0,1.0))) 
    #addfaceset!(grid1, "bottom", (x)->x[2]==0.0)
    #addfaceset!(grid1, "top", (x)->x[2]==1.0)
    addcellset!(grid2, "topcells", (x)->true)

    #f(x) = 0.02*(x-1.5)*(x-11.5)
    #for nodeid in getnodeset(grid2, "bottom")
    #    x = grid2.nodes[nodeid].x
    #    dy = f(x[1])
    #    new_node = Node(Vec{2,Float64}((x[1], x[2]+dy)))
    #    grid2.nodes[nodeid] = new_node
    #end
    grid = gridmerge(grid1, grid2)
    addcellset!(grid, "all", collect(1:getncells(grid)))
    return grid
end

data = ProblemData(
    dim = 2,
    tend = 1.0
)

data.grid = generate_my_grid2()


material = LinearElastic(
    E = 1e5,
    ν = 0.3
)

con1 = Dirichlet(
    set = getfaceset(data.grid, "left2"),
    func = (x,t) -> (0.0, 0.0),
    field = :u,
    dofs = [1,2]
)
push!(data.dirichlet, con1)

con1 = Dirichlet(
    set = getfaceset(data.grid, "right1"),
    func = (x,t) -> (0.0, -1.0),
    field = :u,
    dofs = [1,2]
)
push!(data.dirichlet, con1)

part = Part{2,Float64}(
    element  = Five.SolidElement{2,1,RefCube,Float64}(
        thickness = 1.0, 
        qr_order = 2,
        celltype = Quadrilateral,
        dimstate = PlaneStrain()
    ),
    material = material,
    cellset = sort(collect(getcellset(data.grid, "all")))
)
push!(data.parts, part)

data.output[] = Output(
    interval = data.tend/10,
    runname = "contact",
    savepath = "."
)

masterface = collect( getfaceset(data.grid, "top2")    )
slavenodes = collect( getnodeset(data.grid, "bottom1") )
slavenodes = Five.nodeset_to_vertexset(data.grid, slavenodes)

data.contact = Five.FeSurface{2,Float64}(masterface, slavenodes )

solver = ExplicitSolver(
    Δt0 = 1e-4,
    Δt_max = 0.1,
)

state, data = build_problem(data)

output = solvethis(solver, state, data);