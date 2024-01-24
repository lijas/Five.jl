

function generate_enf_grid(nelx, nely, L, h, a0, SolidCellType, CohesiveCellType)

    grid1 = generate_grid(SolidCellType,(nelx,nely),Vec((0.0,0.0)),Vec((L,h)))
    grid2 = generate_grid(SolidCellType,(nelx,nely),Vec((0.0,h)),Vec((L,h*2)))
    grid = Five.gridmerge(grid1,grid2)

    addvertexset!(grid, "mid", (x)-> x[1] ≈ L/2 && x[2] ≈ h*2)
    @assert(length(getvertexset(grid, "mid")) == 2)
    addvertexset!(grid, "botleft", (x)-> x[1] ≈ 0.0 && x[2] ≈ 0.0)
    addvertexset!(grid, "botright", (x)-> x[1] ≈ L && x[2] ≈ 0.0)

    construct_interfacer_cells!(grid, "top1", "bottom2", CohesiveCellType)

    solid_cells = collect(1:nelx*nely*2)
    cz_cells = collect((1:nelx) .+ 2*nelx*nely)

    addcellset!(grid, "solid_cells", solid_cells)
    addcellset!(grid, "cz_cells", cz_cells)

    addcellset!(grid, "temp_precraced", (x)-> x[1]>L-a0)
    precracked_cells = setdiff(getcellset(grid, "temp_precraced"), solid_cells)

    addcellset!(grid, "precracked", precracked_cells)

    return grid
end

function construct_interfacer_cells!(grid, setname1::String, setname2::String, CohesiveCellType)


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
    for (i,topface_index) in enumerate(grid1_top_faceset)
        botface_index = grid2_bottom_faceset[i]

        topcell = grid.cells[topface_index[1]]
        botcell = grid.cells[botface_index[1]]
        #topface = Ferrite.faces(grid.cells[topface_index[1]])[topface_index[2]]
        #botface = Ferrite.faces(grid.cells[botface_index[1]])[botface_index[2]]
        cz_nodes = Int[]
        for j in [4, 3]
            push!(cz_nodes, topcell.nodes[j])
        end
        for j in [1, 2]
            push!(cz_nodes, botcell.nodes[j])
        end

        #if CohesiveCellType === Five.CZQuadraticQuadrilateral
        #    push!(cz_nodes, topcell.nodes[7])
        #    push!(cz_nodes, botcell.nodes[5])
        #end

        #cell_nodes = [topface[2], topface[1], botface[1], botface[2]]
        new_cell = CohesiveCellType(Tuple(cz_nodes))
        push!(grid.cells, new_cell)
    end

end