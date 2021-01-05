
@testset "Material 2d" begin

    #Check if Material2D produces same result as 2d-version of MaterialElastic

    #
    #Plane strain
    #
    mat2d = Five.MatLinearElastic(E = 10.0, nu = 0.3, plane_stress = false)
    mat3d = Five.MatLinearElastic(E = 10.0, nu = 0.3)
    state2d = Five.getmaterialstate(mat2d)

    mat2d_ = Material2D(mat3d, Five.PLANE_STRAIN)

    ε = rand(SymmetricTensor{2,2})

    s_, ds_, _ = Five.constitutive_driver(mat2d_, ε, state2d)
    s, ds, _ =  Five.constitutive_driver(mat2d, ε, state2d)
    @test all(s .≈ s_)
    @test all(ds .≈ ds_)

    #
    #Plane stress
    #
    mat2d = Five.MatLinearElastic(E = 10.0, nu = 0.3, plane_stress = true)
    mat3d = Five.MatLinearElastic(E = 10.0, nu = 0.3)
    state2d = Five.getmaterialstate(mat2d)

    mat2d_ = Material2D(mat3d, Five.PLANE_STRESS)

    ε = rand(SymmetricTensor{2,2})

    s_, ds_, _ = Five.constitutive_driver(mat2d_, ε, state2d)
    s, ds, _ =  Five.constitutive_driver(mat2d, ε, state2d)

    @test all(s .≈ s_)
    @test all(ds .≈ ds_)

end
