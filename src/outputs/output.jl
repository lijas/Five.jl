export Output, OutputData, VTKOutput
export VTKCellOutput, VTKNodeOutput
export push_output!

abstract type AbstractOutput end
abstract type StateOutput <:AbstractOutput end

"""
    OutputData  


"""
struct OutputData{output <: AbstractOutput} <: AbstractOutput
    type::output
    set#::Ferrite.IndexSets
    data::Vector{Any} #Stores the data for each timestep.
    interval::Float64 
    output_times::Vector{Float64}
end

function OutputData(; type, set, interval)
    return OutputData(type, set, Any[], interval, Float64[])
end

function build_outputdata(output::OutputData{OutputType}, dh) where OutputType
    o = build_outputdata(output.type, output.set, dh)
    return OutputData(o, output.set, output.data, output.interval, output.output_times)
end

#Pass-through function
function collect_output!(output::OutputData, state::StateVariables, globaldata)
    data = collect_output!(output.type, state, output.set, globaldata)
    push!(output.data, data)
end

should_output(o::OutputData, t) =  o.interval <= (t-last(o.output_times))
set_last_output!(o::OutputData, t) = push!(o.output_times, t)

"""
    VTKNodeOutput

Defines the data that should be added to the nodes of the VTK outputfile
"""

struct VTKNodeOutput{output <: AbstractOutput}
    type::output
    func::Function
end

function VTKNodeOutput(; type, func = mean)
    return VTKNodeOutput(type, func)
end

outputname(o::VTKNodeOutput) = outputname(o.type)

"""
VTKCellOutput

Defines the data that should be added to the cells of the VTK outputfile
"""
struct VTKCellOutput{output <: AbstractOutput}
    type::output
    func::Function
end

function VTKCellOutput(; type, func = mean)
    return VTKCellOutput(type, func)
end

outputname(o::VTKCellOutput) = outputname(o.type)

abstract type AbstractVTKOutputType end
struct FerriteVTKOutput <: AbstractVTKOutputType end
struct FiveVTKOutput <: AbstractVTKOutputType end

struct VTKOutput{OT <: AbstractVTKOutputType}
    #data::Union{Any, Vector{Any}}
    pvd::WriteVTK.CollectionFile
    interval::Float64 

    celloutputs::Vector{VTKCellOutput}
    nodeoutputs::Vector{VTKNodeOutput}

    last_output::Base.RefValue{Float64}
end

function VTKOutput{T}(interval::Float64, savepath::String, filename::String) where T
    pvd = paraview_collection(joinpath(savepath, filename))
    return VTKOutput{T}(pvd, interval, VTKCellOutput[], VTKNodeOutput[], Ref(-Inf))
end

mutable struct Output{T}
    savepath::String
    runname::String
    termination::SimulationTermination
    
    vtkoutput::VTKOutput

    export_rawdata::Bool
    export_vtk::Bool

    #List of stuff that should be outputed
    outputdata::Dict{String, AbstractOutput}
end

function Output(; savepath=raw".", runname::String, interval::T, export_rawdata=false, export_vtk=true, vtkoutputtype::Type{<:AbstractVTKOutputType} = FerriteVTKOutput) where {T}

    vtkoutput = VTKOutput{vtkoutputtype}(interval, savepath, runname)

    #Create logging file
    logfile = open(joinpath(savepath, "log.txt"), "w+")
    global LOG = Logging.SimpleLogger(logfile)

    #Create interact.txt file
    try
        file = open(joinpath(savepath, "interact.txt"), "w") 
        close(file)
    catch
        error("Error, could not create file interact.txt")
    end

    return Output{T}(savepath, runname, _NO_TERMINATION, vtkoutput, export_rawdata, export_vtk, Dict{String, OutputData}())
end

function Ferrite.close!(output::Output, dh)
    for (key, outp) in output.outputdata
        #Overwright with new data
        output.outputdata[key] = build_outputdata(outp, dh)
    end
end

set_simulation_termination!(output::Output, status::SimulationTermination) = output.termination = status
get_simulation_termination(output::Output) = output.termination

function _check_output(oh::Output, name::String, output)
    haskey(oh.outputdata, name) && throw(ArgumentError("there already exists a output with the name: $name"))
end

function push_output!(output::Output, name::String, o::AbstractOutput)
    _check_output(output, name, o)
    output.outputdata[name] = o
end

function push_vtkoutput!(output::Output, o::VTKCellOutput)
    push!(output.vtkoutput.celloutputs, o)
end

function push_vtkoutput!(output::Output, o::VTKNodeOutput)
    push!(output.vtkoutput.nodeoutputs, o)
end


function save_outputs(output::Output)
    filename = output.runname * string("_stateoutputs.jld2")
    save(joinpath(output.savepath, filename), output.outputdata)
    vtk_save(output.vtkoutput.pvd)
end

function should_output(output::Union{VTKOutput}, t::T) where T
    return output.interval <= (t-output.last_output[])
end

function vtk_add_state!(output::Output{T}, state::StateVariables, globaldata) where {T}
    if should_output(output.vtkoutput, state.t)
        output.vtkoutput.last_output[] = state.t
        filename =  output.runname * string(state.step)
        _vtk_add_state!(output, state, globaldata, filename=filename)
    end
end


function create_vtk_output(vtkoutput::VTKOutput{FiveVTKOutput}, state::StateVariables, globaldata; filename::String) where {T}
    
    dh    = globaldata.dh
    parts = globaldata.parts
    
    vtmfile = vtk_multiblock(filename)

    #Ouput to vtk_grid
    for (partid, part) in enumerate(parts)
        #Vtk grid
        vtkfile = get_part_vtk_grid(part)
        vtkfile === nothing && continue
        
        multiblock_add_block(vtmfile, vtkfile, "partid$partid")

        #Export Fields, such :u, etc
        for field in get_fields(part)
            data = eval_part_field_data(part, dh, state, field.name)
            vtkfile[string(field.name)] = data
        end

        for nodeoutput in vtkoutput.nodeoutputs
            data = eval_part_node_data(part, nodeoutput, state, globaldata)
            if data !== nothing
                name = outputname(nodeoutput)
                vtk_point_data(vtkfile, data, name)
            end
        end

        for celloutput in vtkoutput.celloutputs
            data = eval_part_cell_data(part, celloutput, state, globaldata)
            if data !== nothing
                name = outputname(celloutput)
                vtk_cell_data(vtkfile, data, name)
            end
        end

    end

    #collection_add_timestep(output.pvd, vtmfile, state.t)
    vtkoutput.pvd[state.t] = vtmfile
end


function export!(output::Output, state::StateVariables, globaldata; force=false)

    #Export vtk/visualization
    if force || should_output(output.vtkoutput, state.t)
        
        output.vtkoutput.last_output[] = state.t
        filename =  output.runname * string(state.step)
        
        if output.export_vtk
            vtkfile = create_vtk_output(output.vtkoutput, state, globaldata, filename=filename)
            vtk_save(vtkfile)
        end

        if output.export_rawdata
            error("Not implemented")
            dumpstate()
        end
    end

    #Export data
    for (name, outp) in output.outputdata
        if force || should_output(outp, state.t)
            set_last_output!(outp, state.t)
            collect_output!(outp, state, globaldata)
        end
    end
end


function create_vtk_output(vtkoutput::VTKOutput{FerriteVTKOutput}, state::StateVariables, globaldata; filename)
    
    (; grid, dh, parts) = globaldata

    vtkgrid = vtk_grid(filename, grid)

    #Export all fields
    vtk_point_data(vtkgrid, dh, state.d)

    #Nodedata
    for nodeoutput in vtkoutput.nodeoutputs
        DT = getdatatype(nodeoutput.type)
        data = [zero(DT)*NaN for _ in 1:getnnodes(grid)]
        for part in parts
            collect_nodedata!(data, part, nodeoutput.type, state, globaldata)
        end
        vtk_point_data(vtkgrid, data, outputname(nodeoutput.type))
    end 
    
    #Celldata
    for celloutput in vtkoutput.celloutputs
        DT = getdatatype(celloutput.type)
        data = [zero(DT)*NaN for _ in 1:getncells(grid)]
        for part in parts
            collect_celldata!(data, part, celloutput.type, state, globaldata)
        end
        vtk_cell_data(vtkgrid, data, outputname(celloutput.type))
    end

    #Partsets
    for (ipart, part) in enumerate(parts)
        data = [NaN for _ in 1:getncells(grid)]
        data[get_cellset(part)] .= 1.0
        vtk_cell_data(vtkgrid, data, "Part $ipart")
    end

    vtkoutput.pvd[state.t] = vtkgrid
    
end


function reshape_vtk_coords(coords::AbstractVector{Vec{dim,T}}) where {dim,T}
    npoints = length(coords)
    out = zeros(T, (dim == 2 ? 3 : dim), npoints)
    out[1:dim, :] .= reshape(reinterpret(T, coords), (dim, npoints))
    return out
end

function outputs!(output::Output{T}, state::StateVariables, globaldata; force = false) where {T}

    for (name, outp) in output.outputdata
        if force || should_output(outp, state.t)
            set_last_output!(outp, state.t)
            collect_output!(outp, state, globaldata)
        end
    end

end


function handle_input_interaction(output)
    #Checkinput
    try
        file = open(joinpath(output.savepath, "interact.txt"), "r+") 
        input = readline(file)
        if input == "abort"
            set_simulation_termination!(output, ABORTED_SIMULATION)
        elseif input == "dump"
            save_outputs(output)
        end

        #Remove everything in file and close
        truncate(file, 0)
        close(file)
    catch e
        println(e)
        error("Could not open file interact.txt")
    end


end



#
struct SubGridGeometry
    grid::Ferrite.AbstractGrid
    cellset::Vector{Int}
    nodemapper::Vector{Int}
    vtkcells::Vector{WriteVTK.MeshCell}
    coords::Matrix{Float64}
end 

function SubGridGeometry(grid::Ferrite.AbstractGrid{dim}, cellset) where dim
    count = 0
    nodemap = Dict{Int,Int}()
    nodemap2 = Int[]
    vtkcells = MeshCell[]
    cellset = collect(cellset)
    subgridnodes = Vec{dim,Float64}[]

    for cellid in cellset
        cell = grid.cells[cellid]
        vtkcelltype = Ferrite.cell_to_vtkcell(typeof(cell))

        nodes = Int[]
        for nodeid in cell.nodes
            nodeid_subgrid = get!(nodemap,nodeid) do 
                push!(subgridnodes, grid.nodes[nodeid].x)
                push!(nodemap2, nodeid)
                count+=1
            end
            push!(nodes, nodeid_subgrid)
        end

        push!(vtkcells, WriteVTK.MeshCell(vtkcelltype, nodes))
    end
    coords = reinterpret(reshape, Float64, subgridnodes)

    return SubGridGeometry(grid, cellset, nodemap2, vtkcells, coords)
end

function WriteVTK.vtk_grid(filename, geometry::SubGridGeometry)
    return vtk_grid(filename, geometry.coords, geometry.vtkcells)
end
