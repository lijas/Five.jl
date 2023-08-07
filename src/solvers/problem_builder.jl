export ProblemData, build_problem

mutable struct ProblemData{dim,T}
    grid::Ferrite.AbstractGrid
    parts::Vector{Five.AbstractPart{dim}}
    dirichlet::Vector{<:Any}
    external_forces::Vector{Five.AbstractExternalForce}
    constraints::Vector{Five.AbstractExternalForce}
    output::Base.RefValue{Output{T}}
    outputdata::Dict{String, Five.AbstractOutput}
    materialstates::Dict{Int, Vector{Any}}
    initial_condition::Vector{Five.InitialCondition}

    t0::T
    tend::T
    adaptive::Bool
end

function ProblemData(; tend::Float64, dim = 3, T = Float64, t0 = 0.0, adaptive = false)
    
    parts = Five.AbstractPart{dim}[]
    dbc   = Any[]
    exfor = Five.AbstractExternalForce[]
    output = Base.RefValue{Output{T}}()
    outputdata = Dict{String, Five.OutputData}()
    cnstr = Five.AbstractExternalForce[]
    states = Dict{Int, Vector{Any}}()
    grid = Grid(Ferrite.AbstractCell[], Node{dim,T}[])
    ic = Five.InitialCondition[]

    return ProblemData{dim,T}(grid, parts, dbc, exfor, cnstr, output, outputdata, states, ic, t0, tend, adaptive)
end

function build_problem(data::ProblemData)
    build_problem((dh, parts, dbc)-> (), data)
end

function build_problem(func!::Function, data::ProblemData{dim,T}) where {dim,T}
    
    _check_input(data)

    #
    dh = MixedDofHandler(data.grid)
    for part in data.parts
        set    = get_cellset(part)
        length(set) == 0 && continue # Some 
        fields = get_fields(part)
        push!(dh, FieldHandler(fields, Set(set))) 
    end 
    close!(dh)
    
    #
    partstates = Vector{AbstractPartState}(undef, getncells(data.grid))
    for part in data.parts
        set    = get_cellset(part)
        states = construct_partstates(part)
        partstates[set] .= states
        #append!(partstates, states)
    end

    for cellid in keys(data.materialstates)
        partstates[cellid].materialstates .= data.materialstates[cellid]
    end
    
    #
    dch = ConstraintHandler(dh)
    for d in data.dirichlet
        add!(dch, d)
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
    ef = ExternalForceHandler{T}(dh)
    for d in data.external_forces
        push!(ef.external_forces, d)
    end
    close!(ef)
    #
    ch = Constraints()
    for d in data.constraints
        push!(ch.external_forces, d)
    end
    close!(ch, dh)

    #output = Output(dh)
    for (name, o) in data.outputdata
        push_output!(data.output[], name, o)
    end
    close!(data.output[], dh) 

    contact = Contact_Node2Segment{dim,T}()# not used

    globaldata = GlobalData(dch, data.grid, dh, ch, ef, contact, data.parts, data.output[], data.t0, data.tend, data.adaptive)

    for part in globaldata.parts
        @show typeof(part)
        init_part!(part, dh)
    end

    func!(dh, data.parts, dch)

    #State
    state = StateVariables(T, ndofs(dh))
    state.t = data.t0
    state.partstates = partstates

    for ic in data.initial_condition
        apply_analytical!(state.d, dh, ic.field, ic.func);
    end

    #System Arrays
    state.system_arrays = SystemArrays(T, ndofs(dh))
    #state.system_arrays.Kⁱ = create_sparsity_pattern(dh, dch)

    return state, globaldata
end

function _check_input(data::ProblemData{dim,T}) where {dim,T}

    all_cellsets = Int[]
    for part in data.parts
        set = get_cellset(part)
        for cellid in set
            if cellid in all_cellsets
                error("$cellid of part $(typeof(part)) in two sets")
            end
        end
        append!(all_cellsets, set)

        !issorted(set) && error("The cellset for the parts must be sorted.")
    end
    @show length(all_cellsets)
    @show  getncells(data.grid)
#    length(all_cellsets) < getncells(data.grid) && error("Not all cells are included in a part.")

end