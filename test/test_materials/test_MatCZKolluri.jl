
function get_MatCZKolluri_loading1(mat)
    jump1 = collect(range(0.0,      stop = mat.δₙ*5, length=100))
    jump2 = collect(range(mat.δₙ*5, stop = 0.0, length=100))
    jump3 = collect(range(0.0,      stop = mat.δₙ*10, length=100))  

    _jump = [jump1..., jump2..., jump3...]
    jump = [Vec{3}((0.0, 0.0, x)) for x in _jump]

    return jump
end

function get_MatCZKolluri_loading2(mat)
    jump1 = collect(range(0.0,      stop = mat.δₜ*5, length=100))
    jump2 = collect(range(mat.δₜ*5, stop = 0.0, length=100))
    jump3 = collect(range(0.0,      stop = mat.δₜ*10, length=100))  

    _jump = [jump1..., jump2..., jump3...]
    jump = [Vec{3}((x/5, x/5, x/5)) for x in _jump]

    return jump
end

function get_MatCZKolluri_loading3(mat)
    jump1 = collect(range(0.0,      stop = mat.δₜ*5, length=100))
    jump2 = collect(range(mat.δₜ*5, stop = 0.0, length=100))
    jump3 = collect(range(0.0,      stop = mat.δₜ*10, length=100))  

    _jump = [jump1..., jump2..., jump3...]
    jump = [Vec{2}((x/2, x/2)) for x in _jump]

    return jump
end

@testset "MatCZKolluri checksum" begin
    material = MatCZKolluri(σₘₐₓ = 10.0, τₘₐₓ = 10.0, Φₙ = 1.0, Φₜ = 1.0)    

    check_checksum(material, get_MatCZKolluri_loading1(material), "MatCZKolluri1")
    check_checksum(material, get_MatCZKolluri_loading2(material), "MatCZKolluri2")
    check_checksum(material, get_MatCZKolluri_loading3(material), "MatCZKolluri3")
end

@testset "MatCZKolluri" begin
    δₙ=0.00001
    G1 = 1.0
    σₘₐₓ = 10.0
    mat = MatCZKolluri(σₘₐₓ = σₘₐₓ, τₘₐₓ = σₘₐₓ, Φₙ = G1, Φₜ = G1)
    state = Five.getmaterialstate(mat, 0.0)

    δvec = []
    append!(δvec, range(0.0, stop=mat.δₙ*3, length=199))
    append!(δvec, range(mat.δₙ*3, stop=0.0, length=100))
    append!(δvec, range(0.0, stop=mat.δₙ*4, length=100))
    
    traction = Float64[]
    for δ in δvec
        τ, dτ, state = Five.constitutive_driver(mat, Vec((0.0, δ, 0.0)), state)
        push!(traction, τ[2])
    end

    τ, dτ, state = Five.constitutive_driver(mat, Vec((0.0, mat.δₙ*4, 0.0)), state)
 
    τ, dτ, state = Five.constitutive_driver(mat, Vec((0.0,0.0, mat.δₙ)), Five.getmaterialstate(mat, 0.0))
    @test σₘₐₓ ≈ τ[3]

    τ, dτ, state = Five.constitutive_driver(mat, Vec((0.0,0.0, 100.0)), Five.getmaterialstate(mat, 0.0))
    @test 0.0 ≈ τ[3]
    @test 1.0 ≈ Five.interface_damage(state, 1)
end