export ProblemData, build_problem

mutable struct ProblemData{dim,T}
    grid::Ferrite.AbstractGrid
    parts::Vector{Five.AbstractPart{dim}}
    dirichlet::Vector{Ferrite.Dirichlet}
    external_forces::Vector{Five.AbstractExternalForce}
    constraints::Vector{Five.AbstractExternalForce}
    output::Base.RefValue{Output{T}}
    outputdata::Dict{String, Five.AbstractOutput}
    materialstates::Dict{Int, Vector{Any}}

    t0::T
    tend::T
    adaptive::Bool
end

function ProblemData(; tend::Float64, dim = 3, T = Float64, t0 = 0.0, adaptive = false)
    
    parts = Five.AbstractPart{dim}[]
    dbc   = Ferrite.Dirichlet[]
    exfor = Five.AbstractExternalForce[]
    output = Base.RefValue{Output{T}}()
    outputdata = Dict{String, Five.OutputData}()
    cnstr = Five.AbstractExternalForce[]
    states = Dict{Int, Vector{Any}}()
    grid = Grid(Ferrite.AbstractCell[], Node{dim,T}[])

    return ProblemData{dim,T}(grid, parts, dbc, exfor, cnstr, output, outputdata, states, t0, tend, adaptive)
end

function build_problem(data::ProblemData)
    build_problem((dh, parts, dbc)-> (), data)
end

function build_problem(func!::Function, data::ProblemData{dim,T}) where {dim,T}
    
    # Create FieldHandlers/cellsets with parts that have the same element-type and fields 
    # This makes it easier to constraints that span multiple parts
    dict = Dict{Pair{Type{<:Ferrite.AbstractCell},Vector{Field}}, Vector{Int}}()
    for part in data.parts
        celltype = Ferrite.getcelltype(data.grid, first(part.cellset))
        fields = get_fields(part)

        #Check if this combination of celltype and fields have already been added
        combo = Pair(celltype, fields)
        if !haskey(dict, combo)
            #Add new cellset
            dict[combo] = part.cellset
        else
            #Combine the cellset
            union!(dict[combo], part.cellset)
        end
    end

    #
    dh = MixedDofHandler(data.grid)
    for (combo, cells) in dict
        fields = last(combo)
        push!(dh, FieldHandler(fields, Set(cells)))
    end 
    close!(dh)
    
    #
    partstates = Vector{AbstractPartState}(undef, getncells(data.grid))
    for part in data.parts
        states = construct_partstates(part)
        partstates[part.cellset] = states
        #append!(partstates, states)
    end

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

    contact = Contact_Node2Segment{dim,T}()# not used

    globaldata = GlobalData(dch, data.grid, dh, ch, ef, contact, data.parts, data.output[], data.t0, data.tend, data.adaptive)

    for part in globaldata.parts
        init_part!(part, dh)
    end

    func!(dh, data.parts, dch)

    #State
    state = StateVariables(T, ndofs(dh))
    state.t = data.t0
    state.partstates = partstates
    state.prev_partstates = deepcopy(partstates)

    #System Arrays
    state.system_arrays = SystemArrays(T, ndofs(dh))
    state.system_arrays.Kⁱ = create_sparsity_pattern(dh)

    return state, globaldata
end