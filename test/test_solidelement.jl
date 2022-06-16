

@testset "test solid element quad" begin
    
    element = SolidElement{2,1,RefCube,Float64}(
        thickness = 1.0, 
        qr_order = 2,
        celltype = Quadrilateral,
        dimstate = PlaneStrain()
    )

    elementstate = Five.EmptyElementState()
    material = MaterialModels.LinearElastic(E=200.0, ν=0.3)
    materialstates = [MaterialModels.initial_material_state(material) for i in 1:4]
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