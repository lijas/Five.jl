

@testset "test solid element calls" begin
    
    element = Five.SolidElement{2,1,RefCube,Float64}(
        thickness = 1.0, 
        qr_order = 2,
        celltype = Quadrilateral,
        dimstate = PlaneStrain()
    )

    elementstate = [Five.EmptyElementState() for i in 1:2^3]
    material = MaterialModels.LinearElastic(E=200.0, ν=0.3)
    materialstates = [MaterialModels.initial_material_state(material) for i in 1:4]
    stresses = [zero(SymmetricTensor{2,3,Float64,6}) for i in 1:4]
    strains = [zero(SymmetricTensor{2,3,Float64,6}) for i in 1:4]

    fe = zeros(Float64, 8)
    ke = zeros(Float64, 8, 8)
    Δue = zeros(Float64, 8)
    due = zeros(Float64, 8)
    Δt = 0.0
    ip = Lagrange{2,RefCube,1}()
    coords = Ferrite.reference_coordinates(ip)

    ue = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0]
    Five.integrate_forcevector!(
        element, 
        elementstate, 
        material, 
        materialstates,
        stresses,
        strains,
        fe, 
        coords, 
        Δue,
        ue,
        due,
        Δt)

    Five.integrate_forcevector_and_stiffnessmatrix!(
        element, 
        elementstate, 
        material, 
        materialstates,
        stresses,
        strains,
        ke, 
        fe, 
        coords, 
        Δue,
        ue,
        due,
        Δt)

    me = zeros(Float64, 8, 8) 
    Five.integrate_massmatrix!(
        element,
        elementstate, 
        material, 
        coords, 
        me, 
        ue, 
        due
    )

    fe2 = zeros(Float64, 8)
    ge = Ref(0.0)
    Five.integrate_dissipation!(
        element,       
        elementstate,  
        material,     
        materialstates, 
        fe2,           
        ge,            
        coords, 
        Δue,           
        ue,            
        due,           
        Δt,            
    ) 

end


@testset "test linear solid element calls" begin
    
    element = Five.LinearSolidElement{2,1,RefCube,Float64}(
        thickness = 1.0, 
        qr_order = 2,
        celltype = Quadrilateral,
        dimstate = PlaneStrain()
    )

    elementstate = [Five.EmptyElementState() for i in 1:2^3]
    material = MaterialModels.LinearElastic(E=200.0, ν=0.3)
    materialstates = [MaterialModels.initial_material_state(material) for i in 1:2^2]
    stresses = [zero(SymmetricTensor{2,3,Float64,6}) for i in 1:2^2]
    strains = [zero(SymmetricTensor{2,3,Float64,6}) for i in 1:2^2]
    fe = zeros(Float64, 8)
    ke = zeros(Float64, 8, 8)
    Δue = zeros(Float64, 8)
    due = zeros(Float64, 8)
    Δt = 0.0
    ip = Lagrange{2,RefCube,1}()
    coords = Ferrite.reference_coordinates(ip)

    ue = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0]
    Five.integrate_forcevector!(
        element, 
        elementstate, 
        material, 
        materialstates,
        stresses,
        strains,
        fe, 
        coords, 
        Δue,
        ue,
        due,
        Δt)

    Five.integrate_forcevector_and_stiffnessmatrix!(
        element, 
        elementstate, 
        material, 
        materialstates,
        stresses,
        strains,
        ke, 
        fe, 
        coords, 
        Δue,
        ue,
        due,
        Δt)

    #=me = zeros(Float64, 8, 8) 
    Five.integrate_massmatrix!(
        element,
        elementstate, 
        material, 
        coords, 
        me, 
        ue, 
        due
    )=#

    #=fe2 = zeros(Float64, 8)
    ge = Ref(0.0)
    Five.integrate_dissipation!(
        element,       
        elementstate,  
        material,     
        materialstates, 
        fe2,           
        ge,            
        coords, 
        Δue,           
        ue,            
        due,           
        Δt,            
    )=# 

end

@testset "test cohesive element calls" begin
    
    material = MatCZBilinear(
        K    = 1.0e5,
        Gᴵ   = (0.5, 0.5, 0.5),
        τᴹᵃˣ = (50.0, 50.0, 50.0),
        η    = 1.0
    ) 

    element = CohesiveElement(
        order = 1,
        thickness = 1.0,
        celltype = CohesiveCell{2,4,2}
    )

    nqp = Ferrite.getnquadpoints(element)
    ndofs = Ferrite.ndofs(element)
    celltype = Ferrite.getcelltype(element)

    elementstate = [Five.EmptyElementState() for i in 1:nqp]
    materialstates = [MaterialModels.initial_material_state(material) for i in 1:nqp]
    stresses = [zero(SymmetricTensor{2,3,Float64,6}) for i in 1:nqp]
    strains = [zero(SymmetricTensor{2,3,Float64,6}) for i in 1:nqp]
    fe = zeros(Float64, ndofs)
    ke = zeros(Float64, ndofs, ndofs)
    Δue = zeros(Float64, ndofs)
    due = zeros(Float64, ndofs)
    Δt = 0.0
    ip = Ferrite.default_interpolation(celltype)
    coords = Ferrite.reference_coordinates(ip)

    ue = zeros(Float64, ndofs)
    Five.integrate_forcevector!(
        element, 
        elementstate, 
        material, 
        materialstates,
        stresses,
        strains,
        fe, 
        coords, 
        Δue,
        ue,
        due,
        Δt)

    Five.integrate_forcevector_and_stiffnessmatrix!(
        element, 
        elementstate, 
        material, 
        materialstates,
        stresses,
        strains,
        ke, 
        fe, 
        coords, 
        Δue,
        ue,
        due,
        Δt)

    me = zeros(Float64, 8, 8) 
    Five.integrate_massmatrix!(
        element,
        elementstate, 
        material, 
        coords, 
        me, 
        ue, 
        due
    )

    fe2 = zeros(Float64, ndofs)
    ge = Ref(0.0)
    Five.integrate_dissipation!(
        element,       
        elementstate,  
        material,     
        materialstates, 
        fe2,           
        ge,            
        coords, 
        Δue,           
        ue,            
        due,           
        Δt,            
    ) 

end


@testset "test bar element calls" begin
    
    element = Five.BarElement{2}(
        area = 1.0,
    )

    nqp = Ferrite.getnquadpoints(element)
    ndofs = Ferrite.ndofs(element)
    celltype = Ferrite.getcelltype(element)
    
    material = MaterialModels.LinearElastic(E=200.0, ν=0.3)
    elementstate = [Five.EmptyElementState() for i in 1:nqp]
    materialstates = [MaterialModels.initial_material_state(material) for i in 1:nqp]
    stresses = [zero(SymmetricTensor{2,3,Float64,6}) for i in 1:nqp]
    strains = [zero(SymmetricTensor{2,3,Float64,6}) for i in 1:nqp]
    fe = zeros(Float64, ndofs)
    ke = zeros(Float64, ndofs, ndofs)
    Δue = zeros(Float64, ndofs)
    due = zeros(Float64, ndofs)
    Δt = 0.0
    ip = Ferrite.default_interpolation(celltype)
    coords = Ferrite.reference_coordinates(ip)
    coords = [Vec{2,Float64}((x[1],0.0)) for x in coords]

    ue = zeros(Float64, ndofs)
    Five.integrate_forcevector_and_stiffnessmatrix!(
        element, 
        elementstate, 
        material, 
        materialstates,
        stresses,
        strains,
        ke, 
        fe, 
        coords, 
        Δue,
        ue,
        due,
        Δt
    )

end