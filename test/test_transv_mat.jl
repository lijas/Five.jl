@testset "MatTransvLinearElastic" begin

ν₁₂ = 0.25
ν₃₃ = 0.4
G₁₂ = 5.0e4
E₁₁ = 25.0e5
E₃₃ = 1.0e5

#Create material rotate 180 deg
mat180 = MatTransvLinearElastic(E₁₁, E₃₃, ν₁₂, ν₃₃, G₁₂, deg2rad(0.0))#1π)
ε = symmetric(one(Tensor{2,3,Float64}))
σ, C, newstate =  Five.constitutive_driver(mat180, ε)
@test all(mat180.C ≈ C) #Tensors should be equal if α=180

#Rotate 90 and -90, and they should be thesame
mat1 = MatTransvLinearElastic(E₁₁, E₃₃, ν₁₂, G₁₂, G₁₂, π/2)
mat2 = MatTransvLinearElastic(E₁₁, E₃₃, ν₁₂, G₁₂, G₁₂, -π/2)
ε = symmetric(one(Tensor{2,3,Float64}))
σ, C1, newstate =  Five.constitutive_driver(mat1, ε)
σ, C2, newstate =  Five.constitutive_driver(mat2, ε)
@test all(C1 ≈ C2)

#
α = deg2rad(10)
mate = Five.MatLinearElastic{3,Float64}(1.0, E₁₁,ν₁₂)
σ, Ce, newstate =  Five.constitutive_driver(mate, symmetric(one(Tensor{2,3,Float64})))
_R = [cos(α) -sin(α) 0; 
          sin(α)  cos(α) 0; 
          0            0 1]
          
R = Tensor{2,3}(Tuple(_R))

#C = symmetric(otimesu(R,R) ⊡ Ce ⊡ otimesu(R',R'))
_A2 = zeros(3,3,3,3)
for i in 1:3, j in 1:3, k in 1:3, l in 1:3, m in 1:3, n in 1:3, o in 1:3, p in 1:3
    _A2[i,j,k,l] += R[i,m] * R[j,n] * R[k,o] * R[l,p] * Ce[m,n,o,p]
end

A2 = symmetric(Tensor{4,3}(Tuple(_A2)))
tovoigt(A2)
#z-value should not change of rotating around z-axis
#tovoigt(mat.C)
#tovoigt(C)


end
