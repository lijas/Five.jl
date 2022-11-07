

@testset "Constraints" begin

    data = ProblemData(
        dim = 3,
        tend = 1.0
    )
    
    data.grid = generate_grid(Hexahedron, (2,2,2), Vec((0.0,0.0,0.0)), Vec((2.0,2.0,4.0)))
    addvertexset!(data.grid, "topright", (x) -> x[1] == 2.0 && x[2] == 2.0 && x[3] == 4.0)
    addcellset!(data.grid, "all", (x) -> true)
    
    material = LinearElastic(
        E = 100.0,
        ν = 0.3,
    )
    
    part = Part{3,Float64}(
        element  = Five.LinearSolidElement{3,1,RefCube,Float64}(
            qr_order = 2,
            celltype = Hexahedron,
        ),
        material = material,
        cellset = collect(1:getncells(data.grid))
    )
    push!(data.parts, part)

    #Traction force bottom
    fc = Five.FollowerConstraint(
        field = :u,
        faces = getfaceset(data.grid, "top"),
        mastervertex = getvertexset(data.grid, "topright") |> first,
        component = 1
    )
    push!(data.dirichlet, fc)
    
    data.output[] = Output(
        interval = -1.0,
        runname = "--",
        savepath = "."
        )
        
    state, globaldata = build_problem(data)
    
    #What to test?
    @test 1==1
end
