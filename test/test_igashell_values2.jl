using Five
using IGA

function check()
    #parameters
    dim = 3; T = Float64; dim_p = dim-1; L = 10.0; h = 0.1; b = 1.0;
    G1 = 1;   τ = 1e1;    λ_f = 2*G1/τ;   K = 1e4;   λ_0 = τ/K
    orders = (3,3); r = 3;
    mesh = IGA.generate_nurbsmesh((6,6),orders,(L,b),dim)
    angles = deg2rad.([0,0,0,0,0,0])
    nlayers = length(angles)

    #materials
    layermats = [Five.MatLinearElastic{3,T}(1.0, 200, 0.3) for i in 1:nlayers]
    interfacematerial = Five.MatCohesive{dim,T,1}(λ_0,λ_f,K,τ)
    layer_mats = Five.LayeredMaterial(layermats, angles)

    #layeredata
    data = Five.IGAShellData{dim}(layer_mats, 
                                        interfacematerial, 
                                        (orders...,r), 
                                        mesh.knot_vectors,
                                        h, nlayers, 
                                        4, 4)
        
    #extraction operator
    Ce_mat, _ = IGA.compute_bezier_extraction_operators(data.orders[1:dim_p]..., data.knot_vectors[1:dim_p]...)
    Ce_vec = IGA.bezier_extraction_to_vectors(Ce_mat)


    cellid = 2
    coords = JuAFEM.getcoordinates(mesh, cellid)

    intdata = Five.IGAShellIntegrationData(data, Ce_vec)

    cv = intdata.cell_values_layered

    for cellid in 1:size(mesh.IEN,2)#ncells
        
        
        IGA.set_current_cellid!(intdata, cellid)
        Ce = Five.get_current_extraction_operator(intdata)
        bezier_coords = IGA.compute_bezier_points(Ce, coords)
        reinit!(cv, bezier_coords)
        Five._reinit_base_functions!(cv, Ce)
        
        ue = rand(T, size(cv.U,1))*0.001
        
        qp = 0
        for ilay in 1:nlayers
            for oqp in 1:Five.getnquadpoints_ooplane(cv) ÷ nlayers
                for iqp in 1:Five.getnquadpoints_inplane(cv)
                    qp += 1
                    #dU
                    dUdξ = gradient( (ξ) -> Five._calculate_u(cv.inp_ip, cv.oop_ip, Ce, ue, ξ), cv.qr.points[qp])
                    for d in 1:dim
                        dU1 = Five.function_parent_derivative(cv, qp, ue, d)
                        dU2 = dUdξ[:,d]
                        @show dU1 ≈ dU2
                    end
                    #G
                    for d in 1:dim
                        G = Five._calculate_G(cv.inp_ip, bezier_coords, h, cv.qr.points[qp], d)
                        @show cv.G[qp][d] ≈ G
                    end
                    #g
                    for d in 1:dim
                        g1 = Five._calculate_g(cv.inp_ip, cv.oop_ip, Ce, bezier_coords, h, ue, cv.qr.points[qp], d)
                        g2 = cv.G[qp][d] + Five.function_parent_derivative(cv, qp, ue, d)
                        @show g1 ≈ g2
                    end

                end
            end
        end
    end
end