@testset "MatTransvLinearElastic" begin

ν₁₂ = 0.25
ν₂₃ = 0.4
G₁₂ = 5.0e4
E₁₁ = 25.0e5
E₃₃ = 1.0e5

#Create material rotate 180 deg
mat180 = MatTransvLinearElastic(E1 = E₁₁, E2 = E₃₃, ν_12 = ν₁₂, ν_23 = ν₂₃, G_12 = G₁₂, α = deg2rad(0.0))#1π)
ε = symmetric(one(Tensor{2,3,Float64}))
σ, C, newstate =  Five.constitutive_driver(mat180, ε)
@test all(mat180.C ≈ C) #Tensors should be equal if α=180

#Rotate 90 and -90, and they should be thesame
#mat1 = MatTransvLinearElastic(E1 = E₁₁, E2 = E₃₃, ν_12 = ν₁₂, G_12 = G₁₂, α = deg2rad(π/2))
#mat2 = MatTransvLinearElastic(E1 = E₁₁, E2 = E₃₃, ν_12 = ν₁₂, G_12 = G₁₂, α = deg2rad(-π/2))
#ε = symmetric(one(Tensor{2,3,Float64}))
#σ, C1, newstate =  Five.constitutive_driver(mat1, ε)
#σ, C2, newstate =  Five.constitutive_driver(mat2, ε)
#@test all(C1 ≈ C2)


end
