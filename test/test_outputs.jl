

function build_simple_problem()

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
    
    con1 = Dirichlet(
        set = getfaceset(data.grid, "top"),
        func = (x,t) -> (0.0),
        field = :u,
        dofs = [3]
    )
    push!(data.dirichlet, con1)
    
    con1 = Dirichlet(
        set = getfaceset(data.grid, "back"),
        func = (x,t) -> (0.0),
        field = :u,
        dofs = [1]
    )
    push!(data.dirichlet, con1)
    
    con1 = Dirichlet(
        set = getfaceset(data.grid, "right"),
        func = (x,t) -> (0.0),
        field = :u,
        dofs = [2]
    )
    push!(data.dirichlet, con1)
    
    part = Part{3,Float64}(
        element  = Five.LinearSolidElement{3,1,RefCube,Float64}(
            qr_order = 2,
            celltype = Hexahedron,
        ),
        material = material,
        cellset = collect(1:getncells(data.grid))
    )
    push!(data.parts, part)
    
    data.output[] = Output(
        interval = -0.1,
        runname = "Beamexample",
        savepath = "."
    )
    
    output = OutputData(
        type = DofValueOutput(
            field = :u,
            dofs = [2]
        ),
        interval = 0.1,
        set = getfaceset(data.grid, "bottom")
    )
    
    data.outputdata["reactionforce"] = output
    
    vtkoutput = VTKNodeOutput(
        type = MaterialStateOutput(
            field = :σ
        ),
        func = mean,
    )
    
    Five.push_vtkoutput!(data.output[], vtkoutput)
    
    solver = NewtonSolver(
        Δt0 = 1.0,
        Δt_max = 1.0,
    )
    
    state, globaldata = build_problem(data)

    return state, globaldata
end


@testset "Outputs" begin

    state, globaldata = build_simple_problem()

    output = OutputData(
        type = DofValueOutput(
            field = :u,
            dofs = [2]
        ),
        interval = 0.1,
        set = getvertexset(globaldata.grid, "topright")
    )

    output = Five.build_outputdata(output, globaldata.dh)

    #translate all nodes
    fill!(state.system_arrays, 0.0)
    fill!(state.d, 1.337)
    fill!(state.system_arrays.fⁱ, 1.337)
    fill!(state.system_arrays.fᵉ, 13.337)
    Five.collect_output!(output, state, globaldata)

    @test output.data[1].displacement[1] == 1.337
    @test output.data[1].fint[1] == 1.337
    @test output.data[1].fext[1] == 13.337

    #Dof value: Faceset

    output = OutputData(
        type = DofValueOutput(
            field = :u,
            dofs = [2]
        ),
        interval = 0.1,
        set = getfaceset(globaldata.grid, "bottom")
    )

    output = Five.build_outputdata(output, globaldata.dh)

    #translate all nodes
    fill!(state.system_arrays, 0.0)
    fill!(state.d, 1.337)
    fill!(state.system_arrays.fⁱ, 1.337)
    fill!(state.system_arrays.fᵉ, 13.337)
    Five.collect_output!(output, state, globaldata)

    @test sum(output.data[1].displacement) ≈ 1.337 * 9
    @test sum(output.data[1].fint) ≈ 1.337 * 9 #9 nodes
    @test sum(output.data[1].fext) ≈ 13.337 * 9 #9 nodes

    #MaterialStateOutput

    output = OutputData(
        type = Five.StressOutput(),
        interval = 0.1,
        set = getcellset(globaldata.grid, "all")
    )

    output = Five.build_outputdata(output, globaldata.dh)

    Five.collect_output!(output, state, globaldata)

    @test length(getcellset(globaldata.grid, "all")) == length(output.data[1])
    @test all(iszero.(output.data[1]))

end
