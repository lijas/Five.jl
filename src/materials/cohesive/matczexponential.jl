"""
 Van den bosch cohesive law, from KIM
"""
# const MatCZVanDenBosch{dim} = MatCZExponentialLaw{dim,true}
# const MatCZOtherName{dim} = MatCZExponentialLaw{dim,false}

struct MatCZKolluri <: AbstractCohesiveMaterial
    #constant, not updated for every element
    σₘₐₓ::Float64
    τₘₐₓ::Float64
    δₙ::Float64
    δₜ::Float64
    Φₙ::Float64
    Φₜ::Float64
end

function MatCZKolluri( ; σₘₐₓ::Float64, τₘₐₓ::Float64, Φₙ::Float64, Φₜ::Float64)
    δₙ = Φₙ/(σₘₐₓ*exp(1))
    δₜ = Φₜ/(τₘₐₓ*exp(1/2)) 
    return MatCZKolluri(σₘₐₓ, τₘₐₓ, δₙ, δₜ, Φₙ, Φₜ)
end

struct MatCZKolluriState <: AbstractMaterialState
    Δ::Tensor{1, 3, Float64, 3}
    T::Tensor{1, 3, Float64, 3}    
    Δ_max::Tensor{1, 3, Float64, 3}
    d::NTuple{2,Float64}
    d_c::NTuple{2,Float64} #Coupling damage factors
end

function getmaterialstate(m::MatCZKolluri, d::Float64=zero(Float64))
    @assert( 0.0 <= d <= 1.0)

    T = Float64
    Δ_max = Vector{T}(undef, 3)

    #Normal direction
    Δ_max[end] = -m.δₙ*log(1-d)

    #Shear directions
    for i in 1:2
        Δ_max[i] = sqrt(-2m.δₙ*log(1-d))
    end

    return MatCZKolluriState(zero(Vec{3,T}), zero(Vec{3,T}), Vec{3,T}(Tuple(Δ_max)), (d,d), (d,d))
end

is_dissipative(::MatCZKolluri) = true

get_material_state_type(::MatCZKolluri) = MatCZKolluriState
n_damage_parameters(m::MatCZKolluri) = 4
function interface_damage(m::MatCZKolluriState, d::Int)
    d == 1 && return m.d[1]
    d == 2 && return m.d[2]
    d == 3 && return m.d_c[1]
    d == 4 && return m.d_c[2]
end
onset_displacement(mat::MatCZKolluri, d::Int)  = (3==d) ? mat.δₙ : mat.δₜ
max_traction_force(mat::MatCZKolluri, d::Int)  = (3==d) ? mat.σₘₐₓ : mat.τₘₐₓ

function _MatCZKolluri_law_with_damage(m::MatCZKolluri, Δ::Vec{3,T}, ms::MatCZKolluriState) where T <: Real
    
    # unpack tangential and normal directions from separation tensor
    Δₜ = Vec{2,T}( (Δ[1], Δ[2]) )
    Δₙ = Δ[3]

    # compute max separations
    Δₜ_max = Vec{2,T}( (max(norm(Δ[1]), ms.Δ_max[1]), max(norm(Δ[2]), ms.Δ_max[2])) )
    Δₙ_max = max(Δ[3], ms.Δ_max[3])
    
    # compute damage variables
    d_n =  (Δₙ_max ==Inf) ? T(1.0) : 1 - exp(-Δₙ_max/m.δₙ)
    d_cn = (Δₙ_max ==Inf) ? T(1.0) : 1 - exp(-Δₙ_max/m.δₙ)*(1 + Δₙ_max/m.δₙ)
    d_ct = (norm(Δₜ_max) ==Inf) ? T(1.0) : 1 - exp(-dot(Δₜ_max,Δₜ_max)/(2m.δₜ^2))
    d_t =  (norm(Δₜ_max) ==Inf) ? T(1.0) : 1 - exp(-dot(Δₜ_max,Δₜ_max)/(2m.δₜ^2))

    #define Heavyside function
    H(Δ) = Δ < 0 ? 0.0 : 1.0

    # compute normal and tangential tractions
    Tₜ = m.Φₜ/m.δₜ^2 * (1-d_t) * (1-d_cn) * Δₜ
    Tₙ  = m.Φₙ/m.δₙ^2 * (1-d_n*H(Δₙ)) * (1-d_ct*H(Δₙ)) * Δₙ
    
    
    if Δₙ_max == Inf || Δₜ_max == Inf
        ΔD = 0.0
    else
        ⁿΔₜ_max = Vec{2,T}( (ms.Δ_max[1], ms.Δ_max[1]) )
        #Dissipaiton from magrnus, (d_cn -> d_tn     d_ct -> d_nt)
        Δd_n = 1/m.δₙ * exp(-Δₙ_max/m.δₙ) * (Δₙ_max - ms.Δ_max[3])
        Δd_cn  = Δₙ_max/m.δₙ^2 * exp(-Δₙ_max/m.δₙ) * (Δₙ_max - ms.Δ_max[3])
        Δd_t  = 1/m.δₜ^2 * exp(-dot(Δₜ_max,Δₜ_max)/(2m.δₜ^2)) * dot(Δₜ_max, (Δₜ_max - ⁿΔₜ_max))
        Δd_ct = Δd_t

        ΔD = 0.5 * (m.Φₙ/m.δₙ^2) * Δₙ^2 * (1 - d_ct) * Δd_n + 
             0.5 * (m.Φₙ/m.δₙ^2) * Δₙ^2 * (1 - d_n)  * Δd_ct + 
             0.5 * (m.Φₜ/m.δₜ^2) * dot(Δₜ,Δₜ) * (1 - d_t) * Δd_cn + 
             0.5 * (m.Φₜ/m.δₜ^2) * dot(Δₜ,Δₜ) * (1 - d_cn)* Δd_t
    end
        

    #return traction tensor
    return Vec{3,T}(((Tₜ[1], Tₜ[2], Tₙ))), d_n, d_t, d_cn, d_ct, Vec{3,T}(((Δₜ_max[1], Δₜ_max[2], Δₙ_max))), ΔD
end


function constitutive_driver(m::MatCZKolluri, J::Vec, ms::MatCZKolluriState)
    J_dual = Tensors._load(J, nothing)
    _T, _d_n, _d_t, _d_cn, _d_ct, _J_max_temp, _ = _MatCZKolluri_law_with_damage(m, J_dual, ms)

    d_n =  Tensors._extract_value(_d_n)
    d_t =  Tensors._extract_value(_d_t)
    d_cn =  Tensors._extract_value(_d_cn)
    d_ct =  Tensors._extract_value(_d_ct)
    J_max_temp = Tensors._extract_value(_J_max_temp)

    T, dTdΔ = Tensors._extract_value(_T), Tensors._extract_gradient(_T, J)

    return T, dTdΔ, MatCZKolluriState(J, T, Vec{3,Float64}(Tuple(J_max_temp)), (d_n,d_t), (d_cn, d_ct))
end


function constitutive_driver(m::MatCZKolluri, _J::Vec{2}, ms::MatCZKolluriState) 

    J = Vec{3,Float64}((_J[1], 0.0, _J[2]))
    _T::Vec{3,Float64}, _dTdΔ::Tensor{2,3,Float64,9}, new_state = constitutive_driver(m, J, ms)

    #Remove third direction
    T = Vec{2,Float64}((_T[1], _T[3]))
    dTdΔ = Tensor{2,2,Float64,4}((_dTdΔ[1,1], _dTdΔ[3,1], _dTdΔ[1,3], _dTdΔ[3,3]))
    
    return T, dTdΔ, new_state
end

function constitutive_driver_dissipation(mp::MatCZKolluri, J::Vec{3}, prev_state::MatCZKolluriState)

    J_dual = Tensors._load(J, nothing)
    _, _, _, _, _, _, _ΔD = _MatCZKolluri_law_with_damage(mp, J_dual, prev_state)

    ΔD, dΔDdJ =  Tensors._extract_value(_ΔD), Tensors._extract_gradient(_ΔD, J)

    return  ΔD, dΔDdJ
end

function constitutive_driver_dissipation(mp::MatCZKolluri, _J::Vec{2}, prev_state::MatCZKolluriState)
    J = Vec{3,Float64}((_J[1], 0.0, _J[2]))
    _ΔD, _dΔDdJ::Vec{3,Float64} = constitutive_driver_dissipation(mp, J, prev_state)

    return _ΔD, Vec{2,Float64}((_dΔDdJ[1], _dΔDdJ[3]))
end