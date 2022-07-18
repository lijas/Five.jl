"""
BarElement

Any dim

"""

struct BarElement{dim} <: AbstractElement
    area::Float64
    fields::Vector{Field}
end

function BarElement{dim}(; area::Float64) where {dim}
    @assert( dim != 1 )
    return BarElement{dim}(area, [Field(:u, Lagrange{1,RefCube,1}(), dim)])
end

#getquadraturerule(e::SolidElement) = QuadratureRule{1,RefCube}(1)
Ferrite.getnquadpoints(::BarElement) = 1
Ferrite.ndofs(::BarElement{dim}) where {dim} = dim*2
Ferrite.nnodes(::BarElement) = 2
Ferrite.getcelltype(::BarElement{dim}) where dim = Cell{dim,2, dim==3 ? 0 : 1}()
has_constant_massmatrix(::BarElement) = true
get_fields(e::BarElement) = return e.fields

function integrate_massmatrix!(element::BarElement{dim}, elstate::AbstractElementState, material::AbstractMaterial, X::Vector{Vec{dim,T}}, me::Matrix, ue::AbstractVector, due::AbstractVector) where {dim,T}
    error("Not implemented")
end

function integrate_forcevector_and_stiffnessmatrix!(element::BarElement{dim}, 
                            elementstate::EmptyElementState, 
                            material::AbstractMaterial, 
                            materialstate::AbstractVector{<:AbstractMaterialState},
                            stresses::Vector{<:SymmetricTensor{2,3,T}},
                            strains::Vector{<:SymmetricTensor{2,3,T}},
                            ke::AbstractMatrix, 
                            fe::Vector{T}, 
                            X, 
                            Δue::Vector,
                            ue::Vector,
                            due::Vector,
                            dt::T) where {dim,T}

    xe = X + reinterpret(Vec{dim,T}, ue)
    A = element.area

    x = reinterpret(Vec{dim,T}, xe)
    
    L = norm(X[1] - X[2])

    #current
    xa = x[1]; xb = x[2]
    l = norm(xa-xb);
    n = (xb-xa)/l;
    λ = l/L;
    ε = log(λ);

    ε_tensor = SymmetricTensor{2,1,T,1}((ε,))
    τ, Et, new_state = material_response(UniaxialStress(), material, ε_tensor, materialstate[1])
    materialstate[1] = new_state
    stresses[1] = dim==3 ? τ        : MaterialModels.increase_dim(τ)
    strains[1]  = dim==3 ? ε_tensor : MaterialModels.increase_dim(ε_tensor)

    dT = exp(-2*ε)*(A/L) * ( (Et[1,1,1,1] - 2*τ[1,1])*(n⊗n) + τ[1,1]*ones(Tensor{2,dim,T}) );
    Tb = A*exp(-ε)*τ[1] * n;

    k = reshape([dT...],dim,dim)
    ke .= vcat(hcat(k,-k),
              hcat(-k,k))
    fe .= [-Tb..., Tb...];

end
