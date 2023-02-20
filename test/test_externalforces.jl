

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
    
    part = Part(
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
            field = :σ,
            datatype = SymmetricTensor{2,3,Float64,6}
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



@testset "External forces" begin

    state, globaldata = build_simple_problem()

    #Traction force bottom
    force = Five.TractionForce(
        set = getfaceset(globaldata.grid, "bottom"),
        traction = (X,t) -> (0.0, 0.0, 1.0), #Force per area
    )

    force = Five.init_external_force!(force, globaldata.dh)
    Five.apply_external_force!(force, state, globaldata)
    @test sum(state.system_arrays.fᵉ) ≈ 4.0 atol = 1e-10

    #Traction force right

    fill!(state.system_arrays, 0.0)
    force = Five.TractionForce(
        set = getfaceset(globaldata.grid, "right"),
        traction = (X,t) -> (1.0, 0.0, 0.0), #Force per area
    )

    force = Five.init_external_force!(force, globaldata.dh)
    Five.apply_external_force!(force, state, globaldata)
    @test sum(state.system_arrays.fᵉ) ≈ 8.0 atol = 1e-10

    #Point force

    fill!(state.system_arrays, 0.0)

    force = PointForce(
        field = :u,
        comps = [3],
        set = getvertexset(globaldata.grid, "topright"),
        func = (X,t) -> (-13.37,)
    )

    force = Five.init_external_force!(force, globaldata.dh)
    Five.apply_external_force!(force, state, globaldata)
    @test sum(state.system_arrays.fᵉ) ≈ -13.37 atol = 1e-10
end
