# # Beam example

using Five
using IGA

function gridmerge(grids::Vararg{IGA.BezierGrid{dim}}) where dim 
    T = Float64
    Cs = getcelltype.(grids)

    #Check if all celltypes are the same
    C = all(first(Cs) .== Cs) ? first(Cs) : Ferrite.AbstractCell
    
    nodes_new = Node{dim,T}[]
    weights_new = T[]
    beos_new = IGA.BezierExtractionOperator{T}[]
    cells_new = C[]

    faceset_new = Dict{String, Set{FaceIndex}}()
    vertexset_new = Dict{String, Set{VertexIndex}}()
    edgeset_new = Dict{String, Set{EdgeIndex}}()
    nodeset_new = Dict{String, Set{Int}}()
    cellset_new = Dict{String, Set{Int}}()

    for (igrid, grid) in enumerate(grids)
        nodeoffset = length(nodes_new)
        celloffset = length(cells_new)

        #For all cells in grid2, increment the nodeids
        for cell in grid.cells
            N = length(cell.nodes)
            M = nfaces(cell)
            CellType = typeof(cell)
            offset_nodes = [n + nodeoffset for n in cell.nodes]
            offset_cell  = CellType(NTuple{N,Int}(tuple(offset_nodes...)));
            push!(cells_new, offset_cell)
        end
        append!(nodes_new, grid.nodes)
        append!(weights_new, grid.weights)
        append!(beos_new, grid.beo)

        #Update sets
        function updata_set(indexset, new_indexset, INDEX)
            grid2_new_facesets = Dict{String, Set{INDEX}}()
            for (facesetname, faceset) in getproperty(grid, indexset)
                tmp_faceset = Set{INDEX}()
                for (cellid, faceidx) in faceset
                    push!(tmp_faceset, INDEX(cellid + celloffset, faceidx))
                end
                new_indexset[facesetname * string(igrid)] =  tmp_faceset
            end
        end

        updata_set(:facesets, faceset_new, FaceIndex)
        updata_set(:vertexsets, vertexset_new, VertexIndex)
        updata_set(:edgesets, edgeset_new, EdgeIndex)

        for (nodesetname, nodeset) in grid.nodesets
            tmp_nodeset = Set{Int}()
            for nodeid in nodeset
                push!(tmp_nodeset, nodeid + nodeoffset)
            end
            nodeset_new[nodesetname * string(igrid)] =  tmp_nodeset
        end

        for (setname, cellset) in grid.cellsets
            tmp_nodeset = Set{Int}()
            for cellid in cellset
                push!(tmp_nodeset, cellid + celloffset)
            end
            cellset_new[setname * string(igrid)] =  tmp_nodeset
        end

    end

    return BezierGrid(cells_new, nodes_new, weights_new, beos_new; 
            cellsets=cellset_new, nodesets=nodeset_new, facesets=faceset_new, edgesets=edgeset_new, vertexsets=vertexset_new)

end

function generate_my_grid()
    nurbs1 = generate_nurbs_patch(:rectangle, (5,5), (2,2); cornerpos=(0.0,0.0), size = (10.0,10.0)) 
    grid1 = BezierGrid(nurbs1)
    addfaceset!(grid1, "bottom", (x)->x[2]==0.0)
    addfaceset!(grid1, "top", (x)->x[2]==10.0)

    nurbs2 = generate_nurbs_patch(:rectangle, (5,5), (2,2); cornerpos=(1.5,10.0), size = (11.5,10.0)) 
    grid2 = BezierGrid(nurbs2)
    addfaceset!(grid2, "bottom", (x)->x[2]==10.0)
    addfaceset!(grid2, "top", (x)->x[2]==20.0)

    grid = gridmerge(grid1,grid2)
    return grid
end

function generate_my_grid2()
    nurbs1 = generate_nurbs_patch(:rectangle, (5,5), (2,2); cornerpos=(1.5,3.5), size = (10.0,10.0)) 
    grid1 = BezierGrid(nurbs1)
    addfaceset!(grid1, "bottom", (x)->x[2]≈ 3.5)
    addfaceset!(grid1, "top", (x)->x[2]≈13.5)
    addnodeset!(grid1, "bottom", (x)->x[2]≈3.5)
    addcellset!(grid1, "stamp", (x)->true)

    nurbs2 = generate_nurbs_patch(:rectangle, (10,4), (2,2); cornerpos=(0.0,0.0), size = (30.0,3.0)) 
    grid2 = BezierGrid(nurbs2)
    addfaceset!(grid2, "bottom", (x)->x[2]==0.0)
    addfaceset!(grid2, "top", (x)->x[2]==3.0)
    addcellset!(grid2, "block", (x)->true)
    

    f(x) = 0.02*(x-1.5)*(x-11.5)
    for nodeid in getnodeset(grid1, "bottom")
        x = grid1.nodes[nodeid].x
        dy = f(x[1])
        new_node = Node(Vec{2,Float64}((x[1], x[2]+dy)))
        grid1.nodes[nodeid] = new_node
    end

    grid = gridmerge(grid1,grid2)
    return grid
end

data = ProblemData(
    dim = 2,
    tend = 2.0
)

data.grid = generate_my_grid2()

material = LinearElastic(
    E = 1e5,
    ν = 0.3
)

con1 = Dirichlet(
    set = getfaceset(data.grid, "bottom2"),
    func = (x,t) -> (0.0, -1.0),
    field = :u,
    dofs = [1,2]
)
push!(data.dirichlet, con1)

con1 = Dirichlet(
    set = getfaceset(data.grid, "top1"),
    func = (x,t) -> (0.0, 0.0),
    field = :u,
    dofs = [1,2]
)
push!(data.dirichlet, con1)

part = Five.IGAPart{2,Float64}(
    element  = Five.SolidElement{2,1,RefCube,Float64}(
        thickness = 1.0, 
        qr_order = 2,
        celltype = getcelltype(data.grid),
        dimstate = PlaneStrain()
    ),
    material = material,
    cellset = sort(collect(1:getncells(data.grid)))
)
push!(data.parts, part)

data.output[] = Output(
    interval = data.tend/10,
    runname = "contactiga",
    savepath = "."
)

#masterface = collect( getfaceset(data.grid, "top2")    )
#slavenodes = sort(collect( getnodeset(data.grid, "bottom1") ))
#slavenodes = Five.nodeset_to_vertexset(data.grid, slavenodes)

#data.contact = Five.FeSurface{2,Float64}(slavenodes, masterface)

solver = ExplicitSolver(
    Δt0 = 1e-4,
    Δt_max = 0.1,
)

state, data = build_problem(data)

output = solvethis(solver, state, data);