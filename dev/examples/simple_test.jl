using Five

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

@testset "External forces" begin

    #Traction force bottom

    force = Five.TractionForce(
        field = :u,
        set = getfaceset(data.grid, "bottom"),
        traction = (X,t) -> (0.0, 0.0, 1.0), #Force per area
        celltype = Hexahedron
    )

    force = Five.init_external_force!(force, globaldata.dh)
    Five.apply_external_force!(force, state, globaldata)
    @test sum(state.system_arrays.fᵉ) ≈ 4.0 atol = 1e-10

    #Traction force right

    fill!(state.system_arrays, 0.0)
    force = Five.TractionForce(
        field = :u,
        set = getfaceset(data.grid, "right"),
        traction = (X,t) -> (1.0, 0.0, 0.0), #Force per area
        celltype = Hexahedron
    )

    force = Five.init_external_force!(force, globaldata.dh)
    Five.apply_external_force!(force, state, globaldata)
    @test sum(state.system_arrays.fᵉ) ≈ 8.0 atol = 1e-10

    #Point force

    fill!(state.system_arrays, 0.0)

    force = PointForce(
        field = :u,
        comps = [3],
        set = getvertexset(data.grid, "topright"),
        func = (X,t) -> (-13.37,)
    )

    force = Five.init_external_force!(force, globaldata.dh)
    Five.apply_external_force!(force, state, globaldata)
    @test sum(state.system_arrays.fᵉ) ≈ -13.37 atol = 1e-10
end

@testset "Outputs" begin

    #Dof value: Vertexset

    output = OutputData(
        type = DofValueOutput(
            field = :u,
            dofs = [2]
        ),
        interval = 0.1,
        set = getvertexset(data.grid, "topright")
    )

    output = Five.build_outputdata(output, globaldata.dh)

    #translate all nodes
    fill!(state.system_arrays, 0.0)
    fill!(state.d, 1.337)
    fill!(state.system_arrays.fⁱ, 1.337)
    fill!(state.system_arrays.fᵉ, 13.337)
    Five.collect_output!(output, state, globaldata)

    @test output.data[1].displacement == 1.337
    @test output.data[1].fint == 1.337
    @test output.data[1].fext == 13.337

    #Dof value: Faceset

    output = OutputData(
        type = DofValueOutput(
            field = :u,
            dofs = [2]
        ),
        interval = 0.1,
        set = getfaceset(data.grid, "bottom")
    )

    output = Five.build_outputdata(output, globaldata.dh)

    #translate all nodes
    fill!(state.system_arrays, 0.0)
    fill!(state.d, 1.337)
    fill!(state.system_arrays.fⁱ, 1.337)
    fill!(state.system_arrays.fᵉ, 13.337)
    Five.collect_output!(output, state, globaldata)

    @test output.data[1].displacement ≈ 1.337
    @test output.data[1].fint ≈ 1.337 * 9 #9 nodes
    @test output.data[1].fext ≈ 13.337 * 9 #9 nodes

    #MaterialStateOutput

    output = OutputData(
        type = Five.StressOutput(),
        interval = 0.1,
        set = getcellset(data.grid, "all")
    )

    output = Five.build_outputdata(output, globaldata.dh)

    Five.collect_output!(output, state, globaldata)

    @test length(getcellset(data.grid, "all")) == length(output.data[1])
    @test all(iszero.(output.data[1]))


end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

