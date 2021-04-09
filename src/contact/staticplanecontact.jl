
struct ContactState
    ξ::Float64
    t_T::Vec{2,T}
end

struct StaticPlaneContact{dim} <: AbstractContactSearchAlgorithm
	
	slaves::Vector{NodeContactEntity}
    
    contact_states::Vector{ContactState}

    n::Vec{dim,Float64}
    a::Vec{dim,Float64}
    xm::Vec{dim,Float64}
    μ::Float64
    c_T::Float64
    penalty::Float64

end

function StaticPlaneContact{dim}(;
    slaves::AbstractVector{NodeContactEntity},
    n::Vec{dim,Float64},
    a::Vec{dim,Float64},
    xm::Vec{dim,Float64},
    μ::Float64,
    c_T::Float64,
    penalty::Float64) where dim


    return StaticPlaneContact{dim}(collect(slaves), ContactState[], n, a, xm, μ, c_t, penalty)

end

function search1!(contact::StaticPlaneContact{dim}, x::AbstractVector) where dim

    contactoutputs = ContactOutput[]

    for slave in contact.slaves 
        contactoutput = search_contact(slave, contact.master, x)

        if contactoutput !== nothing
            push!(contactoutputs, contactoutput)
        end
    end

    return contactoutputs  
end

function update_contact!(contact::StaticPlaneContact{dim}, x::AbstractVector, it::Int = 0) where dim

end

function handle_contact!(contact::StaticPlaneContact{2}, contact_treatment::PenaltyBasedContactWithoutFriction, contactoutputs, f, K, x)
    
    penalty = contact_treatment.penalty
    
    for co in contactoutputs
        slaveid = co.slavenode
        state = contact.contact_states[co.slavenode]
        
        residual(contact, xs, state)

        xs_dual = Tensors._load(xs, nothing)
        _r, _t_t, _ξ = residual(xs_dual)

        fe, ke = Tensors.extract_value(_r), Tensors.extract_gradient(_r)
        t_T = Tensors.extract_value(_t_T)
        ξ   = Tensors.extract_value_ξ

        contact.contact_states[slaveid] = ContactState(ξ, t_T)
    end

end

function residual(contact, xs, state)

    l = norm(contact.a)
    n = contact.n
    xm = contact.xm
    a = contact.a
    ā = a / l

    ξ = (1/l) * (xs - xm) ⋅ ā
    x1 = xm + ξ*ā
    gns = (xs - x1)*n

    p_N = contact.penalty * gns

    Δg_T = (ξ - state.ξ) * ā

    δg_n = n
    δg_T = a
    
    #Friction
    t_tr = state.t_T + contact.c_T*Δg_T
    ϕ = norm(t_tr)   - contact.μ*p_N

    if ϕ <= 0
        t_T = t_tr
    else
        t_T = mat.μ * abs(p_N) * t_tr/norm(t_tr)
    end

    t_T = t_T ⋅ ā

    r = p_N*δg_n + t_T*δg_T

    return r, t_t, ξ
end


function search_contact(slave::NodeContactEntity, master::StaticPlaneContactEntity, x)
    dim = 2
    T = Float64

    xs = Vec{dim,T}(Tuple(x[slave.dofs]))

    _penatration = dot(master.normal, xs-master.origin)

    if _penatration < 0.0
        return ContactOutput(slave, master, -1.0, _penatration, master.normal, -1.0, zero(Vec{dim,T}))
    else
        return nothing
    end

end

