
function read_grid()

    filepath = joinpath(@__DIR__, "phase_grid.mphtxt")
    file = open(filepath, "r")

    nodes = Node{2,Float64}[]
    cells = Quadrilateral[]
    read_cells = false
    counter = 1
    for line in eachline(file)
        counter += 1
        if line == "#Elements"
            read_cells = true
            continue
        end

        if read_cells
            n1,n2,n3,n4 = parse.(Int, String.(split(strip(line))))

            push!(cells, Quadrilateral((n1,n2,n4,n3) .+ 1))
        else

            x,y = parse.(Float64, String.(split(strip(line))))

            node = (x+0.5, y+0.5) |> Vec{2,Float64} |> Node
            push!(nodes, node)
        end

    end

    grid = Grid(cells, nodes)
    addfaceset!(grid, "bottom", (x)-> x[2] ≈ 0.0)
    return grid

end

