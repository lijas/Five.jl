using Five


NELX = 51
NELY = 51

L = 1.0
b = 0.1
a0 = L/2
E = 210.0e3
nu = 0.3
elsize = L/NELY

data = ProblemData(
    dim = 2,
    tend = 1000.0,
)

#grid
data.grid = generate_grid(Quadrilateral, (NELX,NELY), Vec((0.0,0.0)),Vec((L,L)))
addvertexset!(data.grid, "top", x -> x[2] ≈ L)
addcellset!(data.grid, "precracked", x -> x[1] < a0 && (L/2 - elsize*1.5 < x[2] < L/2 + elsize*1.5))

material = 
MatLinearElastic(
    E = 135.0e3,
    nu = 0.18,
    plane_stress = false
)

#
element = PhaseFieldElement{2,1,RefCube,Float64}(
    thickness = 1.0,     
    Gc = 2.7e-3,
    lc = elsize*2,
    μ = E / (2(1+nu)),
    λ = (E*nu) / ((1+nu) * (1 - 2nu)),
    qr_order = 2, 
    celltype = Quadrilateral
)

#
part = Part{2,Float64}(
    element = element,
    material = material,
    cellset  = 1:(NELX*NELY) |> collect
)
push!(data.parts, part)


#Change initial states
for cellid in getcellset(data.grid, "precracked")
    data.elementstates[cellid] = Five.PhaseFieldElementState(element, 1000.0)
end

#
dbc1 = JuAFEM.Dirichlet(
    field = :u,
    set = getfaceset(data.grid, "bottom"),
    func = (x,t)->[0.0, 0.0],
    dofs =  [1,2]
)
push!(data.dirichlet, dbc1)


#
force = PointForce(
    field = :u,
    comps = [2],
    set = getvertexset(data.grid, "top"),
    func = (X,t) -> [1.0]
)
push!(data.external_forces, force)


#
data.output[] = Output(
    interval = -1.0,
    runname = "phase_field",
    savepath = "."
)

#
output = OutputData(
    type = DofValueOutput(
        field = :u,
        dofs = [1,2]
    ),
    interval = -1.0,
    set = getfaceset(data.grid, "top")
)
data.outputdata["reactionforce"] = output

state, globaldata = build_problem(data)

solver = NewtonSolver(
    Δt0 = 0.01,
    Δt_max = 100.0,
)

solver = LocalDissipationSolver(
    Δλ0          = 0.001,
    Δλ_max       = 200.0,
    Δλ_min       = 1e-11,
    ΔL0          = 2.5,
    ΔL_min       = 1e-7,
    ΔL_max       = 10.0,
    sw2d         = 1e-5,
    sw2i         = -Inf,
    optitr       = 5,
    maxitr       = 12,
    maxitr_first_step = 50,
    maxsteps     = 300,
    λ_max        = +2000.0,
    λ_min        = -2000.0,
    tol          = 1e-4,
    max_residual = 1e5,
    finish_criterion = Five.finish_criterion1
)

output = solvethis(solver, state, globaldata)

d = [output.outputdata["reactionforce"].data[i].displacement for i in 1:length(output.outputdata["reactionforce"].data)]
f = [output.outputdata["reactionforce"].data[i].fint for i in 1:length(output.outputdata["reactionforce"].data)]

fig = plot(xlabel = "Displacement", ylabel = "Force")
plot!(fig, d, f, label = "Elastic energy", mark = :o)