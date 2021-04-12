export ProblemData, build_problem

mutable struct ProblemData{dim,T}
    grid::JuAFEM.AbstractGrid
    parts::Vector{Five.AbstractPart{dim}}
    dirichlet::Vector{JuAFEM.Dirichlet}
    external_forces::Vector{Five.AbstractExternalForce}
    constraints::Vector{Five.AbstractExternalForce}
    output::Base.RefValue{Output{T}}
    outputdata::Dict{String, Five.AbstractOutput}
    materialstates::Dict{Int, Vector{Any}}
    contact::Base.RefValue{StaticPlaneContact{dim}}

    t0::T
    tend::T
    adaptive::Bool
end

function ProblemData(; tend::Float64, dim = 3, T = Float64, t0 = 0.0, adaptive = false)
    
    parts = Five.AbstractPart{dim}[]
    dbc   = JuAFEM.Dirichlet[]
    exfor = Five.AbstractExternalForce[]
    output = Base.RefValue{Output{T}}()
    outputdata = Dict{String, Five.OutputData}()
    cnstr = Five.AbstractExternalForce[]
    states = Dict{Int, Vector{Any}}()
    grid = Grid(JuAFEM.AbstractCell[], Node{dim,T}[])
    contact = Base.RefValue{StaticPlaneContact{dim}}()

    return ProblemData{dim,T}(grid, parts, dbc, exfor, cnstr, output, outputdata, states, contact, t0, tend, adaptive)
end

function build_problem(data::ProblemData)
    build_problem((dh, parts, dbc)-> (), data)
end

function build_problem(func!::Function, data::ProblemData{dim,T}) where {dim,T}
    partstates = AbstractPartState[]
    
    dh = MixedDofHandler(data.grid)
    for part in data.parts
        #Add fields to dofhandler
        fields = get_fields(part)
        cells = get_cellset(part)
        push!(dh, FieldHandler(fields, Set(cells)))

        #partstates
        states = construct_partstates(part)
        append!(partstates, states)
    end
    close!(dh)

    for cellid in keys(data.materialstates)
        partstates[cellid].materialstates .= data.materialstates[cellid]
    end
    
    #
    dch = ConstraintHandler(dh)
    for d in data.dirichlet
        fh = getfieldhandler(dh, first(d.faces)[1])
        add!(dch, fh, d)
    end
    close!(dch)
    update!(dch, 0.0)

    #
    #=vbc = ConstraintHandler(dh)
    for d in velocity_constraints
        fh = getfieldhandler(dh, first(d.faces)[1])
        add!(vbc, fh, d)
    end
    close!(vbc)
    update!(vbc, 0.0)=#

    #
    ef = ExternalForceHandler{T}()
    for d in data.external_forces
        push!(ef.external_forces, d)
    end
    close!(ef, dh)
    #
    ch = Constraints()
    for d in data.constraints
        push!(ch.external_forces, d)
    end
    close!(ch, dh)

    for (name, o) in data.outputdata
        push_output!(data.output[], name, o)
    end
    close!(data.output[], dh)
    
    func!(dh, data.parts, dch)

    contact = data.contact[]
    init_contact!(contact, dh)

    x0 = get_x0(dh)

    globaldata = GlobalData(dch, data.grid, dh, ch, ef, contact, data.parts, data.output[], data.t0, data.tend, x0, data.adaptive)

    #State
    state = StateVariables(T, ndofs(dh))
    state.t = data.t0
    state.partstates = partstates
    state.prev_partstates = deepcopy(partstates)

    #System Arrays
    state.system_arrays = SystemArrays(T, ndofs(dh))
    state.system_arrays.Kⁱ = create_sparsity_pattern(dh)

    for part in globaldata.parts
        init_part!(part, dh)
    end

    return state, globaldata
end