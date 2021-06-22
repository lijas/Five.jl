export ImplicitSolver



@with_kw struct ImplicitSolver{dim,T}
    dbc::ConstraintHandler
    vbc::ConstraintHandler

    grid::Grid
    constraints::Constraints
    efh::ExternalForceHandler{dim,T}
    dh::DofHandler
    contact_handler::Vector{ContactHandler}

    parts::Vector{Part}
    materials::Vector{AbstractMaterial}

    cellstates::Vector{AbstractElementState}
    materialstates::Vector{Vector{AbstractMaterialState}}

    endtime::T
    dt::T
    plotinterval::T
    
    v0::Vector{T} #not used yet
    d0::Vector{T} #not used yet

    savefile::String
end


function solvethis(solver::ImplicitSolver{dim,T}) where {dim,T}
	
	#Extract 
    starttime = time()
    Float = T #or Float32

    #Extract
    @unpack_ImplicitSolver solver
	
	#Define variables
    penelty = 20;
	
	#Variables
    beta = 1/4;
    gamma = 1/2;
    NEWTON_TOL = 0.00001;
    MAX_NEWTON_ITR = 10
    newton_itr = 0;
    energy_output_interval = 0.1
    plotinterval_timer = plotinterval;
    energy_output_timer = energy_output_interval
	output = Output{Float}()
	#Some usefull variables
    _ndofs = ndofs(dh)
	#X = getgridcoords(dh);
    update!(dbc,0.0)
    nsteps = convert(Int, round(endtime/dt));
	gravityacc = zeros(T,dim)
    gravityacc[dim] = 0.0#-19.81

    #
    M = create_sparsity_pattern(dh)
    Kint = create_sparsity_pattern(dh)
    C = sparse(zeros(Float, _ndofs, _ndofs))
    f = zeros(Float, _ndofs)
    f_ext_n = zeros(Float, _ndofs)
    f_int_n = zeros(Float, _ndofs)
    f_ext_np1 = zeros(Float, _ndofs)
    f_int_np1 = zeros(Float, _ndofs)
    f_cont = zeros(Float, _ndofs)
    f_constant = zeros(Float, _ndofs)

    #Init matricies
    M = assemble_lumped_massmatrix!(dh, parts, materials, elementinfos, M)

    assemble_bodyforce!(dh, parts, materials, elementinfos, f_constant, gravityacc)
    f_ext_n = copy(f_constant)
    f_ext_np1 = copy(f_constant)

    #Plotting
    #pvd = paraview_collection("C:/Users/elias/Documents/pvd2/mypvd2")
    pvd = paraview_collection(savefile)
    prog = ProgressMeter.Progress(nsteps, 1)
	#prog = ProgressMeter.ProgressThresh(NEWTON_TOL, "Solving:")
    # pre-allocate time stepping variables
	t_n = 0.0;
    d_n = zeros(_ndofs);
    v_n = copy(v0) #v0
    ff = vec(f_ext_n - C*v_n)
    apply_zero!(M, ff, dbc)
    a_n = M\ff

    v_np1 = copy(v0)
    a_np1 = copy(a_n)
    d_np1 = copy(d_n)

    #Energy 
    W_int_np1 = 0.5*(d_np1 - d_n)'*(f_int_n + f_int_np1)
    W_int_n = 0.5*(d_np1 - d_n)'*(f_int_n + f_int_np1)
    W_ext_np1 = 0.5*(d_np1 - d_n)'*(f_ext_n + f_ext_np1)
    W_ext_n = 0.5*(d_np1 - d_n)'*(f_ext_n + f_ext_np1)
    W_kin_np1 = 0.5*v_np1'*M*v_np1

    #f_int_np1, Kint = assemble_stiffnessmatrix_and_forcevector!(dh, parts, materials, materialstates, cellstates, elementinfos, Kint, f_int_np1, d_np1, v_n + (1-gamma)*dt*a_n)

    @timeit "Simulation" for it in 1:nsteps
        t_np1 = it*dt;

        d_tilde_np1 = d_n + dt*v_n + (dt^2)/2*(1 - 2*beta)*a_n;
        v_tilde_np1 = v_n + (1-gamma)*dt*a_n;

        update!(dbc,t_np1)

        d_np1 = copy(d_tilde_np1);

        inewton = 0;
        while true; inewton +=1;

            #Get force                                                                       
            @timeit "Assembling" f_int_np1, Kint = assemble_stiffnessmatrix_and_forcevector!(dh, parts, materials, materialstates, cellstates, elementinfos, Kint, f_int_np1, d_np1, v_tilde_np1)

            f_ext_force = zeros(T, ndofs(dh))
            @timeit "External forces" apply_external_forces!(efh, f_ext_force, t_np1)

            f_ext_np1 = f_constant + f_ext_force
            f = f_ext_np1 - f_int_np1;

            #Search contact
            #contact_constraints = search_node2segment_contact(contact, vec(X+d_np1));
            @timeit "Apply constraint" apply_constraints!(constraints, f, Kint, d_np1);
            
            #applypenalty!(contact_constraints,f,E,vec(X + d_np1));

            #comute a_np1 and v_np1
            a_np1 = 1/(beta*dt^2)*(d_np1 - d_tilde_np1);
            v_np1 = v_tilde_np1 + gamma*dt*a_np1;

            #Compute residule
            r = M*a_np1 - f;

            normg = norm(r[Ferrite.free_dofs(dbc)])

            #Compute jacobian
            A = 1/(beta*dt^2)*M + Kint;

            #Solve 
            apply_zero!(A,-r, dbc)
            @timeit "Solve" dd = A\-r
            apply_zero!(dd, dbc)

            #updata
            d_np1 = d_np1 + dd;

            #Check convergence
            #apply_zero!(f,dbc);
            if normg < NEWTON_TOL
                break;
            end

            if inewton == MAX_NEWTON_ITR
                vtk_save(pvd);
                error("Reached maximum Newton iterations, aborting")
            end
        end
        
        #Check energy balance
        W_int_np1 += 0.5*(d_np1 - d_n)'*(f_int_n + f_int_np1)
        W_ext_np1 += 0.5*(d_np1 - d_n)'*(f_ext_n + f_ext_np1)
        W_kin_np1 = 0.5*v_np1'*M*v_np1
        
        #Update counter
        t_n = t_np1;
        d_n = d_np1;
        v_n = v_np1;
        a_n = a_np1;
        f_int_n = f_int_np1
        f_ext_n = f_ext_np1

        #Output
        @timeit "Save plotfile" if t_n >= plotinterval_timer
            plotinterval_timer += plotinterval;
            vtkfile = vtk_grid(savefile * string(it), dh)
            #collection_add_timestep(pvd, vtk_grid("C:/Users/elias/Documents/pvd/mypvd" * string(plotinterval_timer), dh, vec(d_n)), t_n)
            #collection_add_timestep(pvd, vtk_point_data(vtkfile, dh, vec(d_n), :u), t_n)
            my_point_data = Five_export_vtk(dh, vtkfile, parts, elementinfos, d_n)

            #=cell_eps_vector = T[]
            for cell_eps in materialstates
                maxeps = 0
                for qp_eps in cell_eps
                    if maxeps < qp_eps.ϵp
                        maxeps = qp_eps.ϵp 
                    end
                end
                push!(cell_eps_vector, maxeps)
            end

            my_cell_data = vtk_cell_data(my_point_data, cell_eps_vector, "eps")
            =#
            collection_add_timestep(pvd, my_point_data, t_n)
        end

        if t_n >= energy_output_timer
            energy_output_timer += energy_output_interval 
            push!(output.energy_kinetic, W_kin_np1)
            push!(output.energy_internal, W_int_np1)
            push!(output.energy_external, W_ext_np1)
        end

        ProgressMeter.next!(prog)
        #print("t" * string(t_np1)*"\n")

    end

    endtime = time()
    totaltime = endtime - starttime

    println("number of cells: $(getncells(grid))")
    println("number of timesteps: $(nsteps)")
    println("total time: $(trunc(totaltime/60)) min and $(trunc(totaltime%60)) sec")
    print_timer(linechars = :ascii)

    vtk_save(pvd);
    return output
end
