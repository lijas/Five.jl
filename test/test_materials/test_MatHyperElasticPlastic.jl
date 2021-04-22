 
function get_MatHyperElasticPlastic_loading()
    
    strain1 = range(0.0,  0.02, length=100)
    strain2 = range(0.02, 0.01, length=100)
    strain3 = range(0.01, 0.03, length=100)

    _C = [strain1..., strain2..., strain3...]
    C = [SymmetricTensor{2,3}((x+1, x/10, 0.0, 1.0, 0.0, 1.0)) for x in _C]

    return C
end

@testset "Test MatHyperElasticPlastic" begin
    
    material = MatHyperElasticPlastic(
        elastic_material = MatNeoHook(
            E = 1.0e5,
            ν = 0.3
        ),
        τ₀ = 400.0,
        H = 1.0e5/20
    ) 
    
    loading = get_MatHyperElasticPlastic_loading()

    check_checksum(material, loading, "MatHyperElasticPlastic1")
end