
@testset "tensor op" begin
∂σ∂ɛ = symmetric(rand(Tensor{4,3}))

a = rand(Vec{3,Float64}); 
b = rand(Vec{3,Float64})
a /= norm(a); 
b /= norm(b);
c = Tensors.cross(a,b)
c /= norm(c)
b = Tensors.cross(c,a)

R = Tensor{2,3}(hcat(a,b,c))

A = otimesu(R,R) ⊡ ∂σ∂ɛ ⊡ otimesu(R',R')
B = otimesu(R',R') ⊡ A ⊡ otimesu(R,R)

_A2 = zeros(3,3,3,3)
for i in 1:3, j in 1:3, k in 1:3, l in 1:3, m in 1:3, n in 1:3, o in 1:3, p in 1:3
    _A2[i,j,k,l] += R[i,m] * R[j,n] * R[k,o] * R[l,p] * ∂σ∂ɛ[m,n,o,p]
end

A2 = Tensor{4,3}(Tuple(_A2))

all(A ≈ A2)
all(∂σ∂ɛ ≈ B)

end
