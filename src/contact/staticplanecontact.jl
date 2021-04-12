
struct ContactState
    ξ::Float64
    t_T::Vec{2,Float64}
end

struct StaticPlaneContact{dim} <: AbstractContactSearchAlgorithm
    
    slaveset::Vector{VertexIndex}
	
    n::Vec{dim,Float64}
    a::Vec{dim,Float64}
    xm::Vec{dim,Float64}
    μ::Float64
    c_T::Float64
    penalty::Float64

    slaves::Vector{NodeContactEntity}
    
    in_contact::Vector{Int}
    contact_states::Dict{Int, ContactState}

end

function StaticPlaneContact{dim}(;
    slaveset::Vector{VertexIndex},
    n::Vec{dim,Float64},
    a::Vec{dim,Float64},
    xm::Vec{dim,Float64},
    μ::Float64,
    c_T::Float64,
    penalty::Float64) where dim
    

    return StaticPlaneContact{dim}(collect(slaveset), n, a, xm, μ, c_T, penalty, NodeContactEntity[], Int[], Dict{Int,ContactState}())
end

function init_contact!(contact::StaticPlaneContact{dim}, dh::MixedDofHandler) where dim

    JuAFEM._check_same_celltype(dh.grid, cellid.(contact.slaveset))

    fh = getfieldhandler(dh, cellid(first(contact.slaveset)))
    
    for vertex in contact.slaveset
        dofs = dofs_on_vertex(dh, fh, vertex, :u, collect(1:dim)) |> Tuple |> Vec{dim,Int}
        x0 = vertexcoords(dh, vertex)
        push!(contact.slaves, NodeContactEntity(dofs, x0) )
    end

end

function contact!(contact::StaticPlaneContact, state::StateVariables, x::AbstractVector)
    contactoutputs = search1!(contact, x)
    handle_contact!(contact, contactoutputs, state, x)
end

function search1!(contact::StaticPlaneContact{dim}, x::AbstractVector) where dim

    contactoutputs = ContactOutput[]
    empty!(contact.in_contact)

    for i in eachindex(contact.slaves)
        slave = contact.slaves[i]
        contactoutput = search_contact(slave, contact, x)
        
        if contactoutput !== nothing
            push!(contact.in_contact, i)
            push!(contactoutputs, contactoutput)
        end
    end

    return contactoutputs  
end


function handle_contact!(contact::StaticPlaneContact{2}, contactoutputs, state, x)
    dim = 2
    T = Float64

    for slaveid in contact.in_contact
        
        slave = contact.slaves[slaveid]
        
        xs = slave.X + Vec{dim,T}(Tuple(x[slave.dofs]))
        xs_dual = Tensors._load(xs, nothing)

        if haskey(contact.contact_states, slaveid)
            cstate = contact.contact_states[slaveid]
            _r, _t_t, _ξ = residual(contact, xs_dual, cstate)
            
        else
            _r, _t_t, _ξ = residual2(contact, xs_dual)
        end

        fe, ke = Tensors._extract_value(_r), Tensors._extract_gradient(_r, xs)
        t_T = Tensors._extract_value(_t_t)
        ξ   = Tensors._extract_value(_ξ)
        contact.contact_states[slaveid] = ContactState(ξ, t_T)

        state.system_arrays.Kᵉ[slave.dofs, slave.dofs] .+= -ke
        state.system_arrays.fᵉ[slave.dofs] .+= -fe
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
    gns = (xs - x1) ⋅ n

    p_N = contact.penalty * gns

    Δg_T = (ξ - state.ξ) * ā
    δg_n = n
    δg_T = a
    
    #Friction
    t_tr = state.t_T + contact.c_T*Δg_T
    ϕ = norm(t_tr)   - contact.μ*abs(p_N)

    if ϕ <= 0
        t_T = t_tr
    else
        
        t_T = contact.μ * abs(p_N) * t_tr/norm(t_tr)
    end

    r = p_N*δg_n + (t_T ⋅ ā)*δg_T

    return r, t_T, ξ
end

function residual2(contact, xs)

    l = norm(contact.a)
    n = contact.n
    xm = contact.xm
    a = contact.a
    ā = a / l

    ξ = (1/l) * (xs - xm) ⋅ ā
    x1 = xm + ξ*ā
    gns = (xs - x1)⋅ n

    p_N = contact.penalty * gns
    δg_n = n

    r = p_N*δg_n 

    return r, Vec{2,Float64}((0.0, 0.0)), ξ
end


function search_contact(slave::NodeContactEntity, master::StaticPlaneContact, x)
    dim = 2
    T = Float64

    xs = slave.X + Vec{dim,T}(Tuple(x[slave.dofs]))
    #@show xs
    _penatration = dot(master.n, xs-master.xm)

    if _penatration < 0.0
        return ContactOutput(slave, slave, Vec((-1.0, -1.0)), _penatration, master.n, -1.0, -1.0, zero(Vec{dim,T}), zero(Vec{dim,T}))
    else
        return nothing
    end

end

