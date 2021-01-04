"""
MatCZBilinear

Standard bi-linear cohesive zone law material model, but without coupling in shear and normal direction.

From: 
Camanho, P., & Davila, C. G. (2002). Mixed-Mode Decohesion Finite Elements in for the Simulation Composite of Delamination Materials.

"""

struct MatCZBilinearSingleMode{dim,T} <: AbstractCohesiveMaterial
    δ₀::T #Value of jump when traction is maximum
    δf::T #Complete seperation
    τ_max::T #Maximum traction value
    K_tension::T
    K_compression::T 
end

struct MatCZBilinearSingleModeState{dim,T} <: AbstractMaterialState
    δmax::Vector{T}
    damage::Vector{T}
    traction::Vec{dim,T}
    jump::Vec{dim,T}
end

function MatCZBilinearSingleMode{dim}(δ₀::T, δf::T, τ_max::T, K_traction::T, K_compression::T = K_traction) where {dim,T}
    return MatCZBilinearSingleMode{dim,T}(δ₀, δf, τ_max, K_traction, K_compression)
end

function MatCZBilinearSingleModeState(mp::MatCZBilinearSingleMode{dim,T}, d::Vector{T}=zeros(T,dim)) where {dim,T}
    t = zero(Vec{dim,T}) #Traction vector
    J = zero(Vec{dim,T})
    δmax = zeros(T,dim)
    for i in 1:dim
        δmax[i] = -mp.δ₀*mp.δf/(d[i]*(mp.δf - mp.δ₀) - mp.δf)
        if δmax[i] <= mp.δ₀
            δmax[i] = 0.0
        end
    end
    return MatCZBilinearSingleModeState{dim,T}(δmax, d, t, J)
end

function MatCZBilinearSingleModeState(mp::MatCZBilinearSingleMode{dim,T}, d::T) where {dim,T}
    return MatCZBilinearSingleModeState(mp, fill(d,dim))
end

function MatCZBilinearSingleModeState{dim,T}(mp::MatCZBilinearSingleMode{dim,T}, d::T) where {dim,T}
    return MatCZBilinearSingleModeState(mp, fill(d,dim))
end

get_material_state_type(::MatCZBilinearSingleMode{dim,T}) where {dim,T} = MatCZBilinearSingleModeState{dim,T}

#
interface_damage(state::MatCZBilinearSingleModeState) = state.damage[i]
initial_stiffness(mat::MatCZBilinearSingleMode) = mat.K

max_traction_force(mat::MatCZBilinearSingleMode, ::Int = 1) = mat.τ_max

onset_displacement(mat::MatCZBilinearSingleMode, ::Int = 1) = mat.δ⁰

function _constitutive_driver(mp::MatCZBilinearSingleMode{dim,T}, J::Vec{dim,T2}, prev_state::MatCZBilinearSingleModeState) where {dim,T,T2}
    
    K = mp.K_tension
    delta_f = mp.δf
    delta_0 = mp.δ₀
    delta_max = [T2(prev_state.δmax[i]) for i in 1:dim]
    delta = J

    d = zeros(T2,dim) #damage vector (from 0 to 1)
    t = zeros(T2,dim) #Traction vector

    delta_max[dim] = max(delta_max[dim], delta[dim])
    for i in 1:(dim-1)
        delta_max[i] = max(delta_max[i], abs(delta[i]))
    end

    #If interpenetration, apply "contact-forcd"
    _normal_flag = 0 
    if delta[dim] < 0
        t[dim] = mp.K_compression*delta[dim]
        #If we entered this case, dont deal with the normal-direction in next loop
        _normal_flag = 1
    end

    for i in 1:(dim-_normal_flag)
        if delta_max[i] <= delta_0
            t[i] = K*delta[i]
        elseif delta_0 < delta_max[i] < delta_f
            d[i] = delta_f*(delta_max[i]-delta_0)/(delta_max[i]*(delta_f-delta_0))
            t[i] = (1-d[i])*K*delta[i]      
        elseif delta_max[i] >= delta_f
            for j in 1:dim
                t[j] = 0.0
                d[j] = 1.0
            end          
            break
        end
    end

    return Vec{dim}(Tuple(t)), delta_max, d
end

function constitutive_driver(mp::MatCZBilinearSingleMode{dim,T}, J, prev_state::MatCZBilinearSingleModeState) where {dim,T}
    _, delta_max, damage = _constitutive_driver(mp, J, prev_state)
    dt::Tensor{2,dim,T,dim^2}, t::Vec{dim,T} = JuAFEM.gradient(J -> _constitutive_driver(mp, J, prev_state)[1], J, :all)

    return t, dt, MatCZBilinearSingleModeState{dim,T}(delta_max,damage,t,J)
end
