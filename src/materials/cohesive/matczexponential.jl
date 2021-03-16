"""
 Van den bosch cohesive law, from KIM
"""
# const MatCZVanDenBosch{dim} = MatCZExponentialLaw{dim,true}
# const MatCZOtherName{dim} = MatCZExponentialLaw{dim,false}

struct MatVanDenBosch <: AbstractCohesiveMaterial
    #constant, not updated for every element
    σₘₐₓ::Float64
    τₘₐₓ::Float64
    δₙ::Float64
    δₜ::Float64
    Φₙ::Float64
    Φₜ::Float64
    with_damage::Bool
end

function MatVanDenBosch( ; σₘₐₓ::Float64, τₘₐₓ::Float64, Φₙ::Float64, Φₜ::Float64, with_damage::Bool = true)
    δₙ = Φₙ/(σₘₐₓ*exp(1))
    δₜ = with_damage ? Φₜ/(τₘₐₓ*exp(1/2)) : Φₜ/(τₘₐₓ*sqrt(0.5 * exp(1))) 
    return MatVanDenBosch(σₘₐₓ, τₘₐₓ, δₙ, δₜ, Φₙ, Φₜ, with_damage)
end

struct MatVanDenBoschState <: AbstractMaterialState
    Δ::Tensor{1, 3, Float64, 3}
    T::Tensor{1, 3, Float64, 3}    
    Δ_max::Tensor{1, 3, Float64, 3}
    d::NTuple{2,Float64}
    d_c::NTuple{2,Float64} #Coupling damage factors
end

function getmaterialstate(m::MatVanDenBosch, d::Float64=zero(Float64))
    @assert( 0.0 <= d <= 1.0)

    T = Float64
    Δ_max = Vector{T}(undef, 3)

    #Normal direction
    Δ_max[end] = -m.δₙ*log(1-d)

    #Shear directions
    for i in 1:2
        Δ_max[i] = sqrt(-2m.δₙ*log(1-d))
    end

    return MatVanDenBoschState(zero(Vec{3,T}), zero(Vec{3,T}), Vec{3,T}(Tuple(Δ_max)), (d,d), (d,d))
end

get_material_state_type(::MatVanDenBosch) = MatVanDenBoschState
interface_damage(m::MatVanDenBoschState, d::Int) = m.d[d]
onset_displacement(mat::MatVanDenBosch, d::Int)  = (3==d) ? mat.δₙ : mat.δₜ
max_traction_force(mat::MatVanDenBosch, d::Int)  = (3==d) ? mat.σₘₐₓ : mat.τₘₐₓ

function _vandenbosch_law_with_damage(Δ::Vec{3,T}, ms::MatVanDenBoschState, m::MatVanDenBosch) where T <: Real

    # unpack tangential and normal directions from separation tensor
    Δₜ = Δ[1:2]
    Δₙ = Δ[3]

    # compute max separations
    Δₜ_max = [max(norm(Δ[1]), ms.Δ_max[1]), max(norm(Δ[2]), ms.Δ_max[2])]
    Δₙ_max = max(Δ[end], ms.Δ_max[end])
    
    # compute damage variables
    d_n =  (Δₙ_max ==Inf) ? T(1.0) : 1 - exp(-Δₙ_max/m.δₙ)
    d_cn = (Δₙ_max ==Inf) ? T(1.0) : 1 - exp(-Δₙ_max/m.δₙ)*(1 + Δₙ_max/m.δₙ)
    d_ct = (norm(Δₜ_max) ==Inf) ? T(1.0) : 1 - exp(-Δₜ_max'*Δₜ_max/2m.δₜ^2)
    d_t =  (norm(Δₜ_max) ==Inf) ? T(1.0) : 1 - exp(-Δₜ_max'*Δₜ_max/2m.δₜ^2)

    #@show d_n d_ct d_t d_cn

    #define Heavyside function
    H(Δ) = Δ < 0 ? 0.0 : 1.0

    # compute normal and tangential tractions
    Tₜ = m.Φₜ/m.δₜ^2 * (1-d_t) * (1-d_cn) * Δₜ
    Tₙ  = m.Φₙ/m.δₙ^2 * (1-d_n*H(Δₙ)) * (1-d_ct*H(Δₙ)) * Δₙ
    
    #Dissipaiton from magrnus, (d_cn -> d_tn     d_ct -> d_nt)
    Δd_cn = d_cn - ms.d_c[1] 
    Δd_ct = d_ct - ms.d_c[2]
    Δd_n  = d_n - ms.d[1] 
    Δd_t  = d_t - ms.d[2] 

    ΔD = 0.5 * (m.Φₙ/m.δₙ^2) * Δₙ^2 * (1 - d_ct) * Δd_n + 
         0.5 * (m.Φₙ/m.δₙ^2) * Δₙ^2 * (1 - d_n)  * Δd_ct + 
         1.0 * (m.Φₜ/m.δₜ^2) * dot(Δₜ,Δₜ) * (1 - d_t)*(1 - d_cn) * Δd_cn + 
         1.0 * (m.Φₜ/m.δₜ^2) * dot(Δₜ,Δₜ) * (1 - d_cn)* Δd_t

    #return traction tensor
    return Vec{3,T}(((Tₜ[1], Tₜ[2], Tₙ))), d_n, d_t, d_cn, d_ct, ΔD
end

function _vandenbosch_law_without_damage(Δ::Vec{3}, m::MatVanDenBosch)
    # Δ follows the convention that the first one or two entries correspond to the
    # tangential direction(s) and th last entry corresponds to the normal direction
    Δₜ = Δ[1:2]
    Δₙ = Δ[3]
    Tₜ(Δₜ, Δₙ) = 2*m.Φₜ / m.δₜ^2 * Δₜ * (1 + Δₙ/m.δₙ) * exp(-Δₙ/m.δₙ) * exp(-Δₜ'*Δₜ/m.δₜ^2)
    Tₙ(Δₜ, Δₙ) = m.Φₙ / m.δₙ^2 * Δₙ * exp(-Δₙ/m.δₙ) * exp(-Δₜ'*Δₜ/m.δₜ^2)
    T = Vec{3,Float64}(vcat(Tₜ(Δₜ, Δₙ), Tₙ(Δₜ, Δₙ)))
    
    return T
end

function constitutive_driver(m::MatVanDenBosch, J::Vec, ms::MatVanDenBoschState)
    if m.with_damage
        return _constitutive_driver_with_damage(m, J, ms)
    else
        return _constitutive_driver_without_damage(m, J, ms)
    end
end

function _constitutive_driver_without_damage(m::MatVanDenBosch, J::Vec{3}, ms::MatVanDenBoschState)
    #Compute tractions and gradient of tractions
    dTdΔ, T = gradient(J -> _vandenbosch_law_without_damage(J, m), J, :all)

    return T, dTdΔ, MatVanDenBoschState(J, T, zero(Vec{3,Float64}), (0.0, 0.0), (0.0, 0.0))
end

function _constitutive_driver_with_damage(m::MatVanDenBosch, J::Vec{3}, ms::MatVanDenBoschState)

    J_dual = Tensors._load(J)
    _T, _d_n, _d_t, _d_cn, _d_ct, _ = _vandenbosch_law_with_damage(J_dual, ms, m)

    d_n =  Tensors._extract_value(_d_n)
    d_t =  Tensors._extract_value(_d_t)
    d_cn =  Tensors._extract_value(_d_cn)
    d_ct =  Tensors._extract_value(_d_ct)

    T, dTdΔ = Tensors._extract_value(_T), Tensors._extract_gradient(_T, J)
    
    #
    J_max_temp = Vector{Float64}(undef, 3)
    Δₜ_max = J_max_temp[1:2] = [(norm(J[i]) > ms.Δ_max[i]) ? norm(J[i]) : ms.Δ_max[i] for i=1:2]
    Δₙ_max = J_max_temp[3] = (J[3] > ms.Δ_max[3]) ? J[3] : ms.Δ_max[3]

    return T, dTdΔ, MatVanDenBoschState(J, T, Vec{3,Float64}(Tuple(J_max_temp)), (d_n,d_t), (d_cn, d_ct))
end

function _constitutive_driver_with_damage(m::MatVanDenBosch, _J::Vec{2}, ms::MatVanDenBoschState) 

    J = Vec{3,Float64}((_J[1], 0.0, _J[2]))
    _T::Vec{3,Float64}, _dTdΔ::Tensor{2,3,Float64,9}, new_state = _constitutive_driver_with_damage(m, J, ms)

    #Remove third direction
    T = Vec{2,Float64}((_T[1], _T[3]))
    dTdΔ = SymmetricTensor{2,2,Float64,3}((_dTdΔ[1,1], _dTdΔ[3,1], _dTdΔ[3,3]))
    
    return T, dTdΔ, new_state
end

function constitutive_driver_dissipation(mp::MatVanDenBosch, J::Vec{3}, prev_state::MatVanDenBoschState)
    @assert(mp.with_damage == true)

    J_dual = Tensors._load(J)
    _, _, _, _, _, _ΔD = _vandenbosch_law_with_damage(J_dual, prev_state, mp)

    ΔD, dΔDdJ =  Tensors._extract_value(_ΔD), Tensors._extract_gradient(_ΔD, J)

    return  ΔD, dΔDdJ
end

function constitutive_driver_dissipation(mp::MatVanDenBosch, _J::Vec{2}, prev_state::MatVanDenBoschState)
    J = Vec{3,Float64}((_J[1], 0.0, _J[2]))
    _ΔD, _dΔDdJ::Vec{3,Float64} = constitutive_driver_dissipation(mp, J, prev_state)

    return _ΔD, Vec{2,Float64}((_dΔDdJ[1], _dΔDdJ[3]))
end