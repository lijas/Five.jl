
#=@testset "cohesive material" begin

    K = 1e3
    G = (1.0,1.0,1.0)
    τ = (10.0,10.0,10.0)
    η = 1.0
    mat = MatCZBilinear(K,G,τ,η)

    state = MatCZBilinearState(mat)

    jump = collect(range(0.0, stop=mat.δᶠ[3], length=100))
    jump2 = collect(range(mat.δᶠ[3]/2, stop=mat.δᶠ[3]/4, length=100))
    jump3 = collect(range(mat.δᶠ[3]/2, stop=mat.δᶠ[3], length=100))
    append!(jump, jump2)
    append!(jump, jump3)
    
    δvec = [Vec((0.0,0.0,x)) for x in jump]
    
    traction = Float64[]
    for δ in δvec
        τ, dτ, state = Five.constitutive_driver(mat, δ, state)
        push!(traction, τ[3])
    end

    τ, dτ, state = Five.constitutive_driver(mat, Vec((0.0,0.0, mat.δ⁰[3])), Five.MatCZBilinearState(mat, 0.0))
    @test τ[3] ≈ τ[3]

    τ, dτ, state = Five.constitutive_driver(mat, Vec((0.0,0.0, mat.δᶠ[3])), Five.MatCZBilinearState(mat, 0.0))
    @test 0.0 ≈ τ[3]
    @test 1.0 ≈ Five.interface_damage(state)
end

@testset "cohesive material2" begin

    G1 = 1.0;   τmax = 1.0e1;    λ_f = 2*G1/τmax;   K = 1.0*1e5;   λ_0 = τmax/K
    mat = Five.MatCohesive{3}(λ_0,λ_f,τmax,K)

    state = Five.MatCohesiveState(mat, 0.0)

    jump = range(0.0, stop=mat.δf, length=100)
    δvec = [Vec((0.0,0.0,x)) for x in jump]
    
    traction = Float64[]
    for δ in δvec
        τ, dτ, state = Five.constitutive_driver(mat, δ, state)
        push!(traction, τ[3])
    end

    #plot(jump, traction, reuse=false)

    τ, dτ, state = Five.constitutive_driver(mat, Vec((0.0,0.0, mat.δ₀)), Five.MatCohesiveState(mat, 0.0))
    @test τmax ≈ τ[3]

    τ, dτ, state = Five.constitutive_driver(mat, Vec((0.0,0.0, mat.δf)), Five.MatCohesiveState(mat, 0.0))
    @test 0.0 ≈ τ[3]
    @test 1.0 ≈ state.damage[3]

end=#
function tesdf()
    δₙ=0.00001
    G1 = 1.0
    σₘₐₓ = 10.0
    #mat = MatCZKolluri{3}(σₘₐₓ = σₘₐₓ, τₘₐₓ = σₘₐₓ, Φₙ = G1, Φₜ = G1, with_damage = true)
    mat = MatCZBilinear(
        K    = 1.0e7,
        Gᴵ   = (0.003, 0.003, 0.003).*1,
        τᴹᵃˣ = (60.0, 60.0, 60.0).*1,
        η    = 1.0
    ) 
    state = Five.getmaterialstate(mat, 0.0)
    @show mat.δ⁰
    δvec = []
    append!(δvec, range(0.0, stop=mat.δ⁰[1]*10, length=199))
    append!(δvec, range(mat.δ⁰[1]*10, stop=0.0, length=100))
    append!(δvec, range(0.0, stop=mat.δ⁰[1]*20, length=100))
    

    traction = Float64[]
    for δ in δvec
        τ, dτ, state = Five.constitutive_driver(mat, Vec((0.0, 0.0, δ)), state)
        push!(traction, τ[3])
    end

    fig = plot(δvec, traction, reuse=false)
    display(fig)

    

end

function testg()
    G1 = 2.0
    σₘₐₓ = 10.0
    mat = MatCZKolluri(σₘₐₓ = σₘₐₓ, τₘₐₓ = σₘₐₓ, Φₙ = G1, Φₜ = G1, with_damage = true)
    state = Five.getmaterialstate(mat, 1.0)

    δvec = []
    append!(δvec, range(0.0, stop=mat.δₙ*20, length=100))
    #append!(δvec, range(mat.δₙ*3, stop=0.0, length=100))
    #append!(δvec, range(0.0, stop=mat.δₙ*20, length=100))
    
    traction = Float64[]
    g = 0.0
    for δ in δvec
        deltag, dg = Five.constitutive_driver_dissipation(mat, Vec((0.0, δ, 0.0)), state)
        τ, dτ, state = Five.constitutive_driver(mat, Vec((0.0, δ, 0.0)), state)
        #@show state
        push!(traction, τ[2])
        g += deltag
    end

    @show g
    return δvec, traction
end


