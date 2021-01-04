export ExplicitSolver, Output, solvethis

@with_kw struct ExplicitSolver{dim,T}
	dbc::ConstraintHandler
    vbc::ConstraintHandler

	grid::Grid{dim}
    constraints::Constraints
    linearconstraints::LinearConstraints
    efh::ExternalForceHandler{dim,T}
	dh::DofHandler
	contact_handler::Vector{ContactHandler}


    parts::Vector{Part}
    materials::Vector{AbstractMaterial}

    cellstates::Vector{AbstractElementState}
    materialstates::Vector{Vector{AbstractMaterialState}}

    t0::T # starttime
	endtime::T
	Δt::T
	plotinterval::T
	
    d0::Vector{T} 
	v0::Vector{T}
    a0::Vector{T}

    runpath::String
end


function solvethis(solver::ExplicitSolver{dim,T}) where {dim,T}

    starttime = time()

    #Check for potential errors is input
    check_explicit_solver_input(solver)

	#Extract
    @unpack_ExplicitSolver solver

    #logging and output
    dispfile = joinpath(runpath, "disp")

    logfile = open(joinpath(runpath, "log.txt"), "w+")
    global LOG = Logging.SimpleLogger(logfile)

	#Define variables
    reset_timer!()
    
    energy_output_interval = T(0.1)
    plotinterval_timer = tᵢ + plotinterval;
    energy_output_timer = energy_output_interval

    x0 = get_x0(dh)
    Δd = zeros(T,ndofs(dh))

    #y_coord_dof = dofsnodes[2:2:end]
    update!(dbc,tᵢ)
    update!(vbc,tᵢ)
    #ncons = nconstraints(constraints);
    nsteps = convert(Int, round((endtime-tᵢ)/Δt))
	#output = Output{T}()
    gravityacc = zeros(T,dim)
    gravityacc[dim] = -0.0981

	#
	M = create_sparsity_pattern(dh)
    C = sparse(zeros(T, ndofs(dh), ndofs(dh)))

    #
    cell_groups = construct_threaded_part_assemblers(dh, parts, elementinfos, cellstates, materials, materialstates)
    #Create matricies
    #M = assemble_lumped_massmatrix!(dh, parts, materials, cellstates, elementinfos, M, dᵢ, vᵢ)
    M = assemble_lumped_massmatrix!(dh, parts, materials, cellstates, elementinfos, M, dᵢ, vᵢ)

    assemble_bodyforce!(dh, parts, materials, cellstates, elementinfos, fᶜᵒⁿˢᵗ, gravityacc)

    fᵉˣᵗᵢ = copy(fᵉˣᵗᵢ₊₁)
    fᵉˣᵗᵢ = copy(fᶜᵒⁿˢᵗ)
    
    #Plotting
    prog = ProgressMeter.Progress(nsteps, 1)
	
    # pre-allocate time stepping variables
    f0 = vec(fᵉˣᵗᵢ - C*vᵢ)
    apply_zero!(M, f0, dbc)
    aᵢ = M\f0
    
    #Energy 
    Wⁱⁿᵗᵢ₊₁ = 0.5*(dᵢ₊₁  - dᵢ)'*(fⁱⁿᵗᵢ + fⁱⁿᵗᵢ₊₁)
    Wⁱⁿᵗᵢ = 0.5*(dᵢ₊₁  - dᵢ)'*(fⁱⁿᵗᵢ + fⁱⁿᵗᵢ₊₁)
    Wᵉˣᵗᵢ₊₁ = 0.5*(dᵢ₊₁  - dᵢ)'*(fᵉˣᵗᵢ + fᵉˣᵗᵢ₊₁)
    Wᵉˣᵗᵢ = 0.5*(dᵢ₊₁  - dᵢ)'*(fᵉˣᵗᵢ + fᵉˣᵗᵢ₊₁)
    Wᵏⁱⁿᵢ₊₁ = 0.5*vᵢ₊₁'*M*vᵢ₊₁ 
    #W_potᵢ₊₁ = 0.5*map()

    #min_timesteps = calculate_all_timesteps(dh, parts, materials, elementinfos, dᵢ, vᵢ)  

    #assembling = Assembling{dim,T,Threads.nthreads()}(dh, grid, parts, elementinfos, materials)
    
    for c in contact_handler
        update_contact!(c.search_algo, x0 + dᵢ, 0)
    end

    #Create plot file
    pvd = paraview_collection(dispfile)
    vtkfile = vtk_grid(dispfile * string(0), dh)
    my_point_data = Five_export_vtk(dh, vtkfile, parts, elementinfos, dᵢ)
    collection_add_timestep(pvd, my_point_data, tᵢ)

    Mcopy = deepcopy(M)
    @timeit "Simulation" for it in 1:nsteps

        #Time update
        tᵢ₊₁ = tᵢ + Δt;
        tᵢ₊₀₅  = 0.5*(tᵢ + tᵢ₊₁);
        Δtᵢ₊₀₅  = tᵢ₊₁ - tᵢ

        #Partial update of vel
        @. vᵢ₊₀₅  = vᵢ + (tᵢ₊₀₅  - tᵢ)*aᵢ 

        #Enforce veloctiy bc
        update!(vbc, tᵢ₊₀₅ )
        apply!(vᵢ₊₀₅ , vbc)
        apply_lc!(vᵢ₊₀₅, linearconstraints)

        #Update nodal disp
        @. Δd = Δtᵢ₊₀₅ * vᵢ₊₀₅ 
        @. dᵢ₊₁  = dᵢ + Δd 

        fill!(fⁱⁿᵗᵢ₊₁, 0)  
        fill!(fᵉˣᵗᵢ₊₁, 0)
        
        #@timeit "Assembling" assemble_forcevector!(assembling, fⁱⁿᵗᵢ₊₁, dᵢ₊₁ , vᵢ₊₀₅)  
        #@timeit "Assembling" assemble_forcevector!(dh, parts, materials, materialstates, cellstates, elementinfos,  fⁱⁿᵗᵢ₊₁, dᵢ₊₁ , vᵢ₊₀₅)  
        @timeit "Assembling" begin
            for (ii, cg) in enumerate(cell_groups)
                @timeit "Assembling $ii" assemble_forcevector!(cg, fⁱⁿᵗᵢ₊₁, Δd, dᵢ₊₁, vᵢ₊₀₅)
            end
        end
        
        f_ext_force = zeros(T, ndofs(dh))
        @timeit "External Forces" apply_external_forces!(efh, f_ext_force, tᵢ₊₁)
        fᵉˣᵗᵢ₊₁ = fᶜᵒⁿˢᵗ + f_ext_force

        @. f = fᵉˣᵗᵢ₊₁ - fⁱⁿᵗᵢ₊₁# + f_cont;

        #Search contact
        @timeit "Contact search" begin
            for c in contact_handler
                contact!(c, f, M, x0 + dᵢ₊₁, it)
            end
        end

        #Constraint
        @timeit "Apply contraints" apply_constraints!(constraints, f, M, dᵢ₊₁ )

        #Compute aᵢ₊₁ 
        update!(dbc,tᵢ₊₁)
        
        #Assabmle mass matrix each iteration due to rigid bodies have non-constant mass matrix
        @timeit "Update massmatrix" update_massmatrix!(dh, parts, materials, cellstates, elementinfos, M, dᵢ₊₁, vᵢ₊₀₅)
    
        M = deepcopy(Mcopy)
        apply_lc!(M, f, linearconstraints)

        apply!(M, f, dbc)
        @timeit "Solve" aᵢ₊₁  = M\f
        apply_lc!(aᵢ₊₁, linearconstraints)
        #rd = [265, 266, 267, 268, 269, 270, 271]
        #ar = (1.0,0.0,0.0,0.0) .+aᵢ₊₁[rd[4:7]]
        #@show ar'*ar

        #Second partial update nodal vel
        @. vᵢ₊₁  = vᵢ₊₀₅  + (tᵢ₊₁ - tᵢ₊₀₅ )*aᵢ₊₁ 
        
        #Check energy balance 
        Wⁱⁿᵗᵢ₊₁ = 0.5*(dᵢ₊₁  - dᵢ)'*(fⁱⁿᵗᵢ + fⁱⁿᵗᵢ₊₁)
        Wᵉˣᵗᵢ₊₁ = 0.5*(dᵢ₊₁  - dᵢ)'*(fᵉˣᵗᵢ + fᵉˣᵗᵢ₊₁)
        Wᵏⁱⁿᵢ₊₁ = 0.5*vᵢ₊₁'*M*vᵢ₊₁ 
        
        #Update counter
        tᵢ = tᵢ₊₁
        dᵢ .= dᵢ₊₁ 
        vᵢ .= vᵢ₊₁ 
        aᵢ  .= aᵢ₊₁ 
        fⁱⁿᵗᵢ .= fⁱⁿᵗᵢ₊₁
        fᵉˣᵗᵢ .= fᵉˣᵗᵢ₊₁

        
        #Output
        @timeit "Save plotfile" if tᵢ >= plotinterval_timer
            plotinterval_timer += plotinterval;
            vtkfile = vtk_grid(dispfile * string(it), dh)
            #collection_add_timestep(pvd, vtk_grid("C:/Users/elias/Documents/pvd/mypvd" * string(plotinterval_timer), dh, vec(dᵢ)), tᵢ)
            #collection_add_timestep(pvd, vtk_point_data(vtkfile, dh, vec(dᵢ), :u), tᵢ)

            ϵᵖ_cells = zeros(T, length(grid.cells))
            #get_plastic_strains!(cell_groups[1], ϵᵖ_cells)

            #my_point_data = 
            Five_export_vtk(dh, vtkfile, parts, elementinfos, dᵢ)
            #vtk_cell_data(vtkfile, ϵᵖ_cells, "plastic_strain")
            collection_add_timestep(pvd, vtkfile, tᵢ)

        end

        if tᵢ >= energy_output_timer
            energy_output_timer += energy_output_interval 
            #push!(output.energy_kinetic, Wᵏⁱⁿᵢ₊₁)
            #push!(output.energy_internal, Wⁱⁿᵗᵢ₊₁)
            #push!(output.energy_external, Wᵉˣᵗᵢ₊₁)
        end
        #push!(tmpp, dᵢ₊₁[729])
        ProgressMeter.next!(prog)

        flush(LOG.stream)
    end
    close(LOG.stream)
    
    endtime = time()
    totaltime = endtime - starttime

    print_timer(linechars = :ascii)

    println("Simulation completed")
    vtk_save(pvd);
    
    println("number of cells: $(getncells(grid))")
    println("number of timesteps: $(nsteps)")
    println("total time: $(trunc(totaltime/60)) min and $(trunc(totaltime%60)) sec")


    cellstates, materialstates = _cellgroups_to_vectors(dh, cell_groups)

    #save laststate in the run folder
    # FileIO.save(joinpath(runpath, "state_"*string(tᵢ)*".jld2") cellstates materialstates dᵢ vᵢ aᵢ tᵢ

    return nothing#output
end

function _cellgroups_to_vectors(dh, cellgroups)

    cellstates = Vector{AbstractElementState}()
    materialstates = Vector{Vector{AbstractMaterialState}}()
    
    resize!(cellstates, getncells(dh.grid))
    resize!(materialstates, getncells(dh.grid))
    
    for cellgroup in cellgroups
        for (i,cellid) in enumerate(cellgroup.cellid_mapper)
            cellstates[cellid] = cellgroup.cellstates[i]
            materialstates[cellid] = cellgroup.materialstates[i]
        end
    end

    return cellstates, materialstates
end