export Output, OutputData, VTKOutput
export VTKCellOutput, VTKNodeOutput
export push_output!

abstract type AbstractOutput end
abstract type StateOutput <:AbstractOutput end

"""
    OutputData  


"""
struct OutputData{output <: AbstractOutput}
    type::output
    set#::JuAFEM.IndexSets
    data::Vector{Any} #Stores the data for each timestep.
    interval::Float64 
    last_output::Base.RefValue{Float64}
end

function OutputData(; type, set, interval)
    return OutputData(type, set, Any[], interval, Ref(-Inf))
end

function build_outputdata(output::OutputData{OutputType}, dh) where OutputType
    o = build_outputdata(output.type, output.set, dh)
    return OutputData(o, output.set, output.data, output.interval, output.last_output)
end

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

struct VTKOutput
    #data::Union{Any, Vector{Any}}
    pvd::WriteVTK.CollectionFile
    interval::Float64 

    celloutputs::Vector{VTKCellOutput}
    nodeoutputs::Vector{VTKNodeOutput}

    last_output::Base.RefValue{Float64}
end

function VTKOutput(interval::Float64, savepath::String, filename::String)
    pvd = paraview_collection(joinpath(savepath, filename))
    return VTKOutput(pvd, interval, VTKCellOutput[], VTKNodeOutput[], Ref(-Inf))
end

"""
    SolverStepStatistic

Holds some data relevant for the solvers.
"""
struct SolverStepStatistic
    Δλ::Float64
    ΔL::Float64
    Δt::Float64
    ntries::Int
    n_newton_itr::Int
    success::Bool
end

mutable struct SolverStatisticsOutput
    #solvertype::Type
    step_stats::Vector{SolverStepStatistic}
    nsteps::Int
    totaltime::Float64
end

SolverStatisticsOutput() = SolverStatisticsOutput(SolverStepStatistic[], 0, -1.0)


mutable struct Output{T}
    savepath::String
    runname::String
    termination::SimulationTermination
    
    vtkoutput::VTKOutput
    solverdata::SolverStatisticsOutput

    #List of stuff that should be outputed
    outputdata::Dict{String, OutputData}
end

function Output(; savepath=raw".", runname::String, interval::T) where {T}

    vtkoutput = VTKOutput(interval, savepath, runname)

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

    return Output{T}(savepath, runname, _NO_TERMINATION, vtkoutput, SolverStatisticsOutput(), Dict{String, OutputData}())
end

function JuAFEM.close!(output::Output, dh::MixedDofHandler)
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

function push_output!(output::Output, name::String, o::OutputData)
    _check_output(output, name, o)
    output.outputdata[name] = o
end

function push_vtkoutput!(output::Output, o::VTKCellOutput)
    push!(output.vtkoutput.celloutputs, o)
end

function push_vtkoutput!(output::Output, o::VTKNodeOutput)
    push!(output.vtkoutput.nodeoutputs, o)
end

function WriteVTK.vtk_save(output::Output)
    vtk_save(output.pvd)
end

function save_outputs(output::Output)
    filename = output.runname * string("_stateoutputs.jld2")
    save(joinpath(output.savepath, filename), output.outputdata)
end

function should_output(output::Union{OutputData, VTKOutput}, t::T) where T
    return output.interval <= (t-output.last_output[])
end

function vtk_add_state!(output::Output{T}, state::StateVariables, globaldata) where {T}
    if should_output(output.vtkoutput, state.t)
        output.vtkoutput.last_output[] = state.t
        filename =  output.runname * string(state.step)
        _vtk_add_state!(output, state, globaldata, outputname=filename)
    end
end

function _vtk_add_state!(output::Output{T}, state::StateVariables, globaldata; outputname::String) where {T}
    
    dh::JuAFEM.AbstractDofHandler = globaldata.dh
    parts = globaldata.parts
    
    vtmfile = vtk_multiblock(joinpath(output.savepath, outputname))
    dim = JuAFEM.getdim(dh)

    #Ouput to vtk_grid
    for (partid, part) in enumerate(parts)
        #Vtk grid
        @timeit "grid" cells, coords = get_vtk_grid(dh, part)
        if length(cells) == 0
            continue
        end
        coords = reshape(reinterpret(T, coords), (dim, length(coords)))
        vtkfile = WriteVTK.vtk_grid(vtmfile, coords, cells)
        
        #Displacements 
        @timeit "disp" node_coords = get_vtk_displacements(dh, part, state)
        vtkfile["u"] = reshape_vtk_coords(node_coords)

        @timeit "nodedata" for nodeoutput in output.vtkoutput.nodeoutputs
            data = get_vtk_nodedata(part, nodeoutput, state, globaldata)
            if data !== nothing
                name = string(typeof(nodeoutput.type)) #string(celloutput.name)
                vtk_point_data(vtkfile, data, name)
            end
        end

        @timeit "celldata" for celloutput in output.vtkoutput.celloutputs
            data = get_vtk_celldata(part, celloutput, state, globaldata)
            if data !== nothing
                name = string(typeof(celloutput.type))
                vtk_cell_data(vtkfile, data, name)
            end
        end

    end
    #collection_add_timestep(output.pvd, vtmfile, state.t)
    output.vtkoutput.pvd[state.t] = vtmfile

end

function reshape_vtk_coords(coords::AbstractVector{Vec{dim,T}}) where {dim,T}
    npoints = length(coords)
    out = zeros(T, (dim == 2 ? 3 : dim), npoints)
    out[1:dim, :] .= reshape(reinterpret(T, coords), (dim, npoints))
    return out
end

function outputs!(output::Output{T}, state::StateVariables, globaldata) where {T}

    for (name, outp) in output.outputdata
        if should_output(outp, state.t)
            outp.last_output[] = state.t
            data = collect_output!(outp.type, state, outp.set, globaldata)
            push!(outp.data, data)
        end
    end

end

function output_solverstat!(output::Output{T}, state::StateVariables, solver::AbstractSolver, globaldata) where {T}
    #Store some solver statistics
    stats = SolverStepStatistic(state.Δλ, state.ΔL, state.Δt, state.step_tries, state.newton_itr, state.converged)
    push!(output.solverdata.step_stats, stats)
    output.solverdata.nsteps += 1
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