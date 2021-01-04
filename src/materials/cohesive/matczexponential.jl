"""
 Van den bosch cohesive law, from KIM
"""
# const MatCZVanDenBosch{dim} = MatCZExponentialLaw{dim,true}
# const MatCZOtherName{dim} = MatCZExponentialLaw{dim,false}

struct MatVanDenBosch{dim} <: AbstractCohesiveMaterial
    #constant, not updated for every element
    σₘₐₓ::Float64
    τₘₐₓ::Float64
    δₙ::Float64
    δₜ::Float64
    Φₙ::Float64
    Φₜ::Float64
    with_damage::Bool
end

function MatVanDenBosch{dim}( ; σₘₐₓ::Float64, τₘₐₓ::Float64, Φₙ::Float64, Φₜ::Float64, with_damage::Bool = true) where {dim}
    δₙ = Φₙ/(σₘₐₓ*exp(1))
    δₜ = with_damage ? Φₜ/(τₘₐₓ*exp(1/2)) : Φₜ/(τₘₐₓ*sqrt(0.5 * exp(1))) 
    return MatVanDenBosch{dim}(σₘₐₓ, τₘₐₓ, δₙ, δₜ, Φₙ, Φₜ, with_damage)
end

struct MatVanDenBoschState{dim} <: AbstractMaterialState
    Δ::Tensor{1, dim, Float64, dim}
    T::Tensor{1, dim, Float64, dim}    
    Δ_max::Tensor{1, dim, Float64, dim}
    d::NTuple{2,Float64}
end

function MatVanDenBoschState(m::MatVanDenBosch{dim}, d::Float64=zero(Float64)) where {dim}
    @assert( 0.0 <= d <= 1.0)

    T = Float64
    Δ_max = Vector{T}(undef, dim)

    #Normal direction
    Δ_max[end] = -m.δₙ*log(1-d)

    #Shear directions
    for i in 1:dim-1
        Δ_max[i] = sqrt(-2m.δₙ*log(1-d))
    end

    #To avoid numerical issues change Δ_max from Inf to something big
    for i in 1:dim
        if Δ_max[i] == Inf
            Δ_max[i] = m.δₙ*100
        end
    end

    return MatVanDenBoschState{dim}(zero(Vec{dim,T}), zero(Vec{dim,T}), Vec{dim,T}(Tuple(Δ_max)), (d,d))
end

get_material_state_type(::MatVanDenBosch{dim}) where {dim} = MatVanDenBoschState
interface_damage(m::MatVanDenBoschState, d::Int) = m.d[d]
onset_displacement(mat::MatVanDenBosch{dim}, d::Int) where dim = (dim==d) ? mat.δₙ : mat.δₜ
max_traction_force(mat::MatVanDenBosch{dim}, d::Int) where dim = (dim==d) ? mat.σₘₐₓ : mat.τₘₐₓ

function _vandenbosch_law_with_damage(Δ::Vec{dim,T}, ms::MatVanDenBoschState, m::MatVanDenBosch) where dim where T <: Real

    # unpack tangential and normal directions from separation tensor
    Δₜ = Δ[1:(dim-1)]
    Δₙ = Δ[dim]

    # compute max separations
    Δₜ_max = [max(norm(Δ[i]), ms.Δ_max[i]) for i=1:(dim-1)]
    Δₙ_max = max(Δ[end], ms.Δ_max[end])
    
    # compute damage variables
    d_n = 1 - exp(-Δₙ_max/m.δₙ)
    d_ct = 1 - exp(-Δₜ_max'*Δₜ_max/2m.δₜ^2)
    d_t = 1 - exp(-Δₜ_max'*Δₜ_max/2m.δₜ^2)
    d_cn = 1 - exp(-Δₙ_max/m.δₙ)*(1 + Δₙ_max/m.δₙ)

    #define Heavyside function
    H(Δ) = Δ < 0 ? 0.0 : 1.0

    # compute normal and tangential tractions
    Tₜ = m.Φₜ/m.δₜ^2 * (1-d_t) * (1-d_cn) * Δₜ
    Tₙ  = m.Φₙ/m.δₙ^2 * (1-d_n*H(Δₙ)) * (1-d_ct*H(Δₙ)) * Δₙ
    
    #return traction tensor
    return Vec{dim,T}(Tuple(vcat(Tₜ, Tₙ)))
end

function _vandenbosch_law_without_damage(Δ::Vec{dim}, m::MatVanDenBosch) where dim
    # Δ follows the convention that the first one or two entries correspond to the
    # tangential direction(s) and th last entry corresponds to the normal direction
    Δₜ = Δ[1:end-1]
    Δₙ = Δ[end]
    Tₜ(Δₜ, Δₙ) = 2*m.Φₜ / m.δₜ^2 * Δₜ * (1 + Δₙ/m.δₙ) * exp(-Δₙ/m.δₙ) * exp(-Δₜ'*Δₜ/m.δₜ^2)
    Tₙ(Δₜ, Δₙ) = m.Φₙ / m.δₙ^2 * Δₙ * exp(-Δₙ/m.δₙ) * exp(-Δₜ'*Δₜ/m.δₜ^2)
    T = Tensor{1,dim}(vcat(Tₜ(Δₜ, Δₙ), Tₙ(Δₜ, Δₙ)))
    return T
end

function constitutive_driver(m::MatVanDenBosch, J::Vec, ms::MatVanDenBoschState)
    if m.with_damage
        return _constitutive_driver_with_damage(m, J, ms)
    else
        return _constitutive_driver_without_damage(m, J, ms)
    end
end

function _constitutive_driver_without_damage(m::MatVanDenBosch, J::Vec{dim}, ms::MatVanDenBoschState{dim}) where dim
    #Compute tractions and gradient of tractions
    dTdΔ, T = gradient(J -> _vandenbosch_law_without_damage(J, m), J, :all)
    d_n, d_t = (0.0, 0.0)

    return T, dTdΔ, MatVanDenBoschState{dim}(J, T, Vec{dim,Float64}(ntuple(i-> 0.0, dim)), (d_n,d_t))
end

function _constitutive_driver_with_damage(m::MatVanDenBosch, J::Vec{dim}, ms::MatVanDenBoschState{dim}) where dim
    
    # update material state
    # TO DO: it is ugly to compute this twice (here and in the ts_law function)
    J_max_temp = Vector{Float64}(undef, dim)
    Δₜ_max = J_max_temp[1:end-1] = [(norm(J[i]) > ms.Δ_max[i]) ? norm(J[i]) : ms.Δ_max[i] for i=1:(dim-1)]
    Δₙ_max = J_max_temp[end] = (J[end] > ms.Δ_max[end]) ? J[end] : ms.Δ_max[end]
    #@show Δₙ_max
    d_n = 1 - exp(-Δₙ_max/m.δₙ)
    d_t = 1 - exp(-Δₜ_max'*Δₜ_max/2m.δₜ^2)

    #compute tractions and gradient of tractions
    dTdΔ, T = gradient(J->_vandenbosch_law_with_damage(J, ms, m), J, :all)
    
    return T, dTdΔ, MatVanDenBoschState{dim}(J, T, Vec{dim,Float64}(Tuple(J_max_temp)), (d_n,d_t))
end
