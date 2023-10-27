
mutable struct ProblemData{dim,T}
    grid::Ferrite.AbstractGrid
    parts::Vector{Five.AbstractPart}
    constraints_ferrite::Vector{<:Any}
    external_forces::Vector{Five.AbstractExternalForce}
    constraints_five::Vector{Five.AbstractExternalForce}
    outputdata::Dict{String, Five.AbstractOutput}
    vtk_output::Vector{Union{VTKNodeOutput, VTKCellOutput}}
    materialstates::Dict{Int, Vector{Any}}
    # initial_condition::Vector{Five.InitialCondition} 

    runname::String
    savepath::String
    vtkoutputtype::AbstractVTKOutputType
    vtk_output_interval::T
    t0::T
    tend::T
    adaptive::Bool
end

function ProblemData(; runname::String, savepath = ".", tend::Float64, dim = 3, T = Float64, t0 = 0.0, adaptive = false, vtk_output_interval=0.0, vtkoutputtype = FiveVTKOutput())
    
    grid = Grid(Ferrite.AbstractCell[], Node{dim,T}[])
    parts = Five.AbstractPart{dim}[]
    dbc   = Any[]
    exfor = Five.AbstractExternalForce[]
    cnstr = Five.AbstractExternalForce[]
    outputdata = Dict{String, Five.OutputData}()
    vtkoutput = Vector{Union{VTKNodeOutput, VTKCellOutput}}()
    states = Dict{Int, Vector{Any}}()
    #ic = Five.InitialCondition[]

    return ProblemData{dim,T}(grid, parts, dbc, exfor, cnstr, outputdata, vtkoutput, states, runname, savepath, vtkoutputtype, vtk_output_interval, t0, tend, adaptive)
end

function build_problem(data::ProblemData)
    build_problem((dh, parts, dbc)-> (), data)
end

function build_problem(func!::Function, data::ProblemData{dim,T}) where {dim,T}
    
    _check_input(data)

    nparts = length(data.parts)

    #
    dh = DofHandler(data.grid)
    for part in data.parts
        set = get_cellset(part)
        sdh = SubDofHandler(dh, Set(set))
        length(set) == 0 && continue 
        for (field_name, field_ip) in get_fields(part)
            add!(sdh, field_name, field_ip)
        end
        init_part!(part, dh)
    end
    close!(dh)
    
    #
    partstates = Vector{AbstractPartState}()
    for (partid, part) in enumerate(data.parts)
        state = construct_partstates(part)
        push!(partstates, state)
    end

    #
    for cellid in keys(data.materialstates)
        found_part = false
        for partid in 1:nparts
            partcells = get_cellset(data.parts[partid])
            if insorted(cellid, partcells)
                found_part = true
                lcellid = searchsortedfirst(partcells, cellid)
                partstates[partid].materialstates[lcellid] .= data.materialstates[cellid]
            end
        end
        @assert found_part
    end
    
    #
    dch = ConstraintHandler(dh)
    for d in data.constraints_ferrite
        add!(dch, d)
    end
    close!(dch)

    #
    ef = ExternalForceHandler{T}(dh)
    for d in data.external_forces
        push!(ef.external_forces, d)
    end
    close!(ef)

    #
    ch = Constraints()
    for d in data.constraints_five
        push!(ch.external_forces, d)
    end
    close!(ch, dh)

    output = Output(; 
        savepath=data.savepath,
        runname = data.runname, 
        interval = data.vtk_output_interval, 
        #export_rawdata=false, 
        #export_vtk=true, 
        vtkoutputtype = FerriteVTKOutput
    )
    for o in data.vtk_output
        push_vtkoutput!(output, o)
    end
    
    for (name, o) in data.outputdata
        push_output!(output, name, o)
    end
    close!(output, dh) 

    contact = Contact_Node2Segment{dim,T}() #not used

    globaldata = GlobalData(data.grid, dh, dch, ch, ef, contact, data.parts, output, data.t0, data.tend, data.adaptive)

    for part in globaldata.parts
        init_part!(part, dh)
    end

    func!(dh, data.parts, dch)

    #State
    state = StateVariables(T, ndofs(dh))
    state.t = data.t0
    state.partstates = partstates

    #
    #for ic in data.initial_condition
    #    apply_analytical!(state.d, dh, ic.field, ic.func);
    #end

    #System Arrays
    state.system_arrays = SystemArrays(T, ndofs(dh))

    return state, globaldata
end

function _check_input(data::ProblemData{dim,T}) where {dim,T}

    all_cellsets = Int[]
    for (partid, part) in enumerate(data.parts)
        set = get_cellset(part)
        for cellid in set
            if cellid in all_cellsets
                error("$cellid is in two parts (one of which is part $(partid), type = $(typeof(part)))")
            end
        end
        append!(all_cellsets, set)

        !issorted(set) && error("The cellset for the part $(partid) is not sorted.")
    end

    #length(all_cellsets) < getncells(data.grid) && error("Not all cells are included in a part.")
end

