
function calculate_element_volume(cell_values::Five.IGAShellValues, nlayers::Int)
    qp = 0
    V = 0
    for qp in 1:Five.getnquadpoints(cell_values)
        V += getdetJdV(cell_values, qp)
    end
    return V
end

function calculate_element_area(cell_values::Five.IGAShellValues, INDEX)
    A = 0
    for qp in 1:Five.getnquadpoints(cell_values)
        A += Five.getdetJdA(cell_values, qp, INDEX)
    end
    return A
end

function get_and_reinit_cv(igashell, grid, cellid)
    cv = Five.build_cellvalue!(igashell, cellid)
    Ce = Five.get_extraction_operator(Five.intdata(igashell), cellid)
    IGA.set_bezier_operator!(cv, Ce)
    
    coords = getcoordinates(grid, cellid)
    bezier_coords = IGA.compute_bezier_points(Ce, coords)
    
    Five.reinit!(cv, bezier_coords)
    Five.build_nurbs_basefunctions!(cv)
    return cv
end

function get_and_reinit_fv(igashell, grid, index)
    cv = Five.build_facevalue!(igashell, index)
    Ce = Five.get_extraction_operator(Five.intdata(igashell), index[1])
    IGA.set_bezier_operator!(cv, Ce)

    coords = getcoordinates(grid, index[1])
    bezier_coords = IGA.compute_bezier_points(Ce, coords)
    
    Five.reinit!(cv, bezier_coords)
    Five.build_nurbs_basefunctions!(cv)
    return cv
end

function get_curved_mesh(cellstate; h, b, R)

    dim = 3; T = Float64
    ν = 0.4; E = 1e5
    visc_para = 0.0
    
    #Mesh
    orders = (3,3); r = 3 
    angles = deg2rad.(T[0,90,0])
    nlayers = length(angles)
    ninterfaces = nlayers-1
    
    nelx = 100; nely = 1
    nurbsmesh = Five.IGA.generate_curved_nurbsmesh((nelx,nely), orders, pi/2, R, b, multiplicity=(1,1))
    grid = Five.IGA.convert_to_grid_representation(nurbsmesh)
    
    cellstates = [cellstate for i in 1:nelx*nely]
    adaptive = false
    linear = true

    interface_damage = [0.0 for _ in 1:ninterfaces, _ in 1:nelx*nely]
    
    #Material
    interfacematerial = Five.MatCohesive{dim,T}(-1.0, -1.0, -1.0, 1.0, 1.0)
    layermats = [MatLinearElastic{dim}(1.0, E, ν) for i in 1:nlayers]
    layer_mats = Five.LayeredMaterial(layermats, angles)
    
    
    #IGAshell
    data = Five.IGAShellData{dim}(layer_mats, interfacematerial, visc_para, (orders...,r), nurbsmesh.knot_vectors, h, dim==2 ? b : 1.0, nlayers, cellstates, interface_damage, adaptive, linear, 4, 3, 4)
    igashell = Five.IGAShell(collect(1:getncells(grid)), reverse(nurbsmesh.IEN, dims=1), data) 

    return grid, igashell
end

@testset "igashellvalues_curved" begin

    # # #
    # TEST IGASHELL VALUES
    # # #

    h = 0.2;    b = .556;     R = 4.22;
    grid, igashell = get_curved_mesh(Five.LAYERED, h=h, b=b, R=R)
    addedgeset!(grid, "left", (x)-> isapprox(x[1], 0.0, atol=1e-5))
    addedgeset!(grid, "right", (x)-> x[3]≈0.0)
    addedgeset!(grid, "front", (x)-> x[2]≈0.0)
    
    rightface = collect(getedgeset(grid, "right"))
    rightedge = [Five.EdgeInterfaceIndex(edgeidx..., 1) for edgeidx in rightface]
    leftface  = collect(getedgeset(grid, "left"))
    leftedge = [Five.EdgeInterfaceIndex(edgeidx..., 2) for edgeidx in leftface]
    frontface  = collect(getedgeset(grid, "front"))
    frontedgetop = [Five.EdgeInterfaceIndex(edgeidx..., 2) for edgeidx in frontface]
    frontedgebot = [Five.EdgeInterfaceIndex(edgeidx..., 1) for edgeidx in frontface]

    #Volume
    V = 0.0
    for cellid in 1:getncells(grid)

        cv = get_and_reinit_cv(igashell, grid, cellid)

        V += calculate_element_volume(cv, Five.nlayers(igashell))
    end

    #Edgelength 1
    L = 0.0
    for edge in rightedge
        cellid, edgeid, interface = edge
        cv = get_and_reinit_fv(igashell, grid, edge)

        L += calculate_element_area(cv, edge)
    end
    @test L ≈ b

    #Edgelength 2
    L = 0.0
    for edge in leftedge
        cellid, edgeid, interface = edge
        cv = get_and_reinit_fv(igashell, grid, edge)

        L += calculate_element_area(cv, edge)
    end
    @test L ≈ b

    #Edgelength 3
    L = 0.0
    for edge in frontedgetop
        cellid, edgeid, interface = edge
        cv = get_and_reinit_fv(igashell, grid, edge)

        L += calculate_element_area(cv, edge)
    end
    @test isapprox(L, (R+h/2)*pi/2, atol=1e-3)
    
    #Edgelength 4

    L = 0.0
    for edge in frontedgebot
        cellid, edgeid, interface = edge
        cv = get_and_reinit_fv(igashell, grid, edge)

        L += calculate_element_area(cv, edge)
    end
    @test isapprox(L, (R-h/2)*pi/2, atol=1e-3)

    #Side area 1
    A = 0.0
    for edge in rightface
        cellid, edgeid = edge
        cv = get_and_reinit_fv(igashell, grid, edge)

        A  += calculate_element_area(cv,edge)
    end
    @test A ≈ h*b

    #Side area 2
    A = 0.0
    for edge in leftface
        cellid, edgeid = edge
        cv = get_and_reinit_fv(igashell, grid, edge)

        A  += calculate_element_area(cv, edge)
    end
    @test A ≈ h*b

    #Side area 3
    A = 0.0
    for edge in frontface
        cellid, edgeid = edge
        cv = get_and_reinit_fv(igashell, grid, edge)

        A  += calculate_element_area(cv, edge)
    end
    @test isapprox(A, R*pi/2 * h, atol=1e-3)

    #Check curvature
    cv = get_and_reinit_cv(igashell, grid, getncells(grid)÷2)
    κ = getindex.(cv.κᵐ,1,1)
    @test all( isapprox.(κ, 1/R, atol=1e-2) )
    
    @test isapprox(V, (R*pi/2 * b * h), atol=1e-1) #approximatlay equal do due mesh discrit



    # # #
    # TEST LOCAL DOF GETTER
    # # #

    @test Five.igashelldofs(igashell, first(frontedgebot)) == [1, 2, 3, 31, 32, 33, 61, 62, 63, 91, 92, 93]

end

@testset "igashell utils" begin

    order = 2;
    ninterfaces = 1
    @test (Five.generate_knot_vector(order, ninterfaces, 0) .≈ Float64[-1,-1,-1, 1,1,1]) |> all
    @test (Five.generate_knot_vector(order, ninterfaces, 3) .≈ Float64[-1,-1,-1, 0,0,0, 1,1,1]) |> all
    @test (Five.generate_knot_vector(order, ninterfaces, 4) .≈ Float64[-1,-1,-1, 0,0,0,0, 1,1,1]) |> all

    ninterfaces = 2
    @test (Five.generate_knot_vector(order, ninterfaces, 1) .≈ Float64[-1,-1,-1, -1/3, 1/3, 1,1,1]) |> all

    ninterfaces = 2
    @test (Five.generate_knot_vector(order, ninterfaces, [1,2]) .≈ Float64[-1,-1,-1, -1/3, 1/3,1/3, 1,1,1]) |> all
    @test (Five.generate_knot_vector(order, ninterfaces, [2,1]) .≈ Float64[-1,-1,-1, -1/3, -1/3,1/3, 1,1,1]) |> all

    ninterfaces = 3
    @test (Five.generate_knot_vector(order, ninterfaces, [1,1,2]) .≈ Float64[-1,-1,-1, -0.5, 0.0, 0.5,0.5, 1,1,1]) |> all

    order = 1; ninterfaces = 3
    @test (Five.generate_knot_vector(Five.LUMPED, order, ninterfaces) .≈ Float64[-1,-1, 1,1]) |> all
    @test (Five.generate_knot_vector(Five.LAYERED, order, ninterfaces) .≈ Float64[-1,-1, -5/10, 0, 5/10, 1,1]) |> all
    @test (Five.generate_knot_vector(Five.WEAK_DISCONTINIUOS_AT_INTERFACE(1), order, ninterfaces) .≈ Float64[-1,-1, -5/10, -5/10, 1,1]) |> all
    @test (Five.generate_knot_vector(Five.WEAK_DISCONTINIUOS_AT_INTERFACE(3), order, ninterfaces) .≈ Float64[-1,-1, 5/10, 5/10, 1,1]) |> all
    @test (Five.generate_knot_vector(Five.STRONG_DISCONTINIUOS_AT_INTERFACE(2), order, ninterfaces) .≈ Float64[-1,-1, -0.5, 0.0,0.0, 0.5, 1,1]) |> all

    ##
    dim_s = 3;
    ooplane_order = 2; 
    nlayers = 2; ninterfaces = nlayers-1;
    @test 9 == Five.ndofs_per_controlpoint(ooplane_order, nlayers, ninterfaces, dim_s, Five.LUMPED)
    @test 15 == Five.ndofs_per_controlpoint(ooplane_order, nlayers, ninterfaces, dim_s, Five.LAYERED)
    @test 18 == Five.ndofs_per_controlpoint(ooplane_order, nlayers, ninterfaces, dim_s, Five.STRONG_DISCONTINIUOS(1))
    @test 18 == Five.ndofs_per_controlpoint(ooplane_order, nlayers, ninterfaces, dim_s, Five.WEAK_DISCONTINIUOS(1))

    ##
    dim_s = 3;
    ooplane_order = 1; 
    nlayers = 4; ninterfaces = nlayers-1;
    @test  6 == Five.ndofs_per_controlpoint(ooplane_order, nlayers, ninterfaces, dim_s, Five.LUMPED)
    @test 15 == Five.ndofs_per_controlpoint(ooplane_order, nlayers, ninterfaces, dim_s, Five.LAYERED)
    @test 18 == Five.ndofs_per_controlpoint(ooplane_order, nlayers, ninterfaces, dim_s, Five.STRONG_DISCONTINIUOS(1))
    @test 18 == Five.ndofs_per_controlpoint(ooplane_order, nlayers, ninterfaces, dim_s, Five.STRONG_DISCONTINIUOS(2))
    @test 21 == Five.ndofs_per_controlpoint(ooplane_order, nlayers, ninterfaces, dim_s, Five.STRONG_DISCONTINIUOS(3))
    @test 18 == Five.ndofs_per_controlpoint(ooplane_order, nlayers, ninterfaces, dim_s, Five.STRONG_DISCONTINIUOS(4))
    @test 12 == Five.ndofs_per_controlpoint(ooplane_order, nlayers, ninterfaces, dim_s, Five.WEAK_DISCONTINIUOS_AT_INTERFACE(3))
    @test 12 == Five.ndofs_per_controlpoint(ooplane_order, nlayers, ninterfaces, dim_s, Five.WEAK_DISCONTINIUOS_AT_INTERFACE(2))
    @test 18 == Five.ndofs_per_controlpoint(ooplane_order, nlayers, ninterfaces, dim_s, Five.WEAK_DISCONTINIUOS(3))

    ## Bezier extraction matrix
    nlayers = 3; order = 1
    knot_lumped = Five.generate_knot_vector(order, nlayers-1, 0)

    #lumped to layered
    new_knots = [-1 + 2i/(nlayers) for i in 1:(nlayers-1)]
    Cmat_lu2la = Five.generate_out_of_plane_extraction_operators(knot_lumped, order, new_knots, fill(order, nlayers-1))

    #lumped to discont
    knot_layered = Five.generate_knot_vector(order, nlayers-1, order)
    Cmat_la2di = Five.generate_out_of_plane_extraction_operators(knot_layered, order, new_knots, [1,0])
    @test all(Cmat_la2di*Cmat_lu2la .≈ [1.0 0.0; 2/3 1/3; 2/3 1/3; 1/3 2/3; 0.0 1.0])

    #Multiplicity vector
    ninterfaces = 4; order = 2
    @test ([3,0,0,0] .== Five.generate_nmultiplicity_vector(Five.WEAK_DISCONTINIUOS(1), ninterfaces, order)) |> all
    @test ([3,2,2,2] .== Five.generate_nmultiplicity_vector(Five.STRONG_DISCONTINIUOS(1), ninterfaces, order)) |> all

    ## Adaptivity stuff
    h = 0.2;    b = .556;     R = 4.22;
    grid, igashell = get_curved_mesh(Five.LAYERED, h=h, b=b, R=R)

    Five.get_upgrade_operator(Five.adapdata(igashell), Five.LUMPED, Five.LAYERED)
    Five.get_upgrade_operator(Five.adapdata(igashell), Five.WEAK_DISCONTINIUOS(1), Five.FULLY_DISCONTINIUOS)
    Five.get_upgrade_operator(Five.adapdata(igashell), Five.STRONG_DISCONTINIUOS(1), Five.FULLY_DISCONTINIUOS)
    Five.get_upgrade_operator(Five.adapdata(igashell), Five.LUMPED, Five.WEAK_DISCONTINIUOS(1))
    #Five.@showm Five.IGA.beo2matrix(C)'

    #Active basefunctions in layer
    order = 2; ilay = 1; ninterfaces=2
    @test (Five.get_active_basefunction_in_layer(1, order, Five.LUMPED) .== 1:order+1) |> all
    @test (Five.get_active_basefunction_in_layer(2, order, Five.LUMPED) .== 1:order+1) |> all

    @test (Five.get_active_basefunction_in_layer(1, order, Five.LAYERED) .== 1:order+1)  |> all
    @test (Five.get_active_basefunction_in_layer(2, order, Five.LAYERED) .== (1:order+1) .+ order) |> all

    @test (Five.get_active_basefunction_in_layer(1, order, Five.FULLY_DISCONTINIUOS) .== (1:order+1)) |> all
    @test (Five.get_active_basefunction_in_layer(2, order, Five.FULLY_DISCONTINIUOS) .== (1:order+1) .+ (order+1)) |> all

    @test (Five.get_active_basefunction_in_layer(2, order, Five.STRONG_DISCONTINIUOS_AT_INTERFACE(2)) .== (1:order+1) .+ (order)) |> all
    @test (Five.get_active_basefunction_in_layer(3, order, Five.STRONG_DISCONTINIUOS_AT_INTERFACE(2)) .== (1:order+1) .+ 5) |> all

    @test (Five.get_active_basefunction_in_layer(1, order, Five.WEAK_DISCONTINIUOS_AT_INTERFACE(1)) .== (1:order+1)) |> all
    @test (Five.get_active_basefunction_in_layer(2, order, Five.WEAK_DISCONTINIUOS_AT_INTERFACE(1)) .== (1:order+1) .+ order.+1) |> all
    @test (Five.get_active_basefunction_in_layer(3, order, Five.WEAK_DISCONTINIUOS_AT_INTERFACE(1)) .== (1:order+1) .+ order.+1) |> all

end