
struct StaticPlaneContact{dim} <: AbstractContactSearchAlgorithm
	
	slaves::Vector{NodeContactEntity}
	master::StaticPlaneContactEntity

end

function search1!(contact::StaticPlaneContact{dim}, x::AbstractVector) where dim

    contactoutputs = ContactOutput[]

    for slave in contact.slaves 
        contactoutput = search_contact(slave, contact.master, x)

        if contactoutput != nothing
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
        N = co.normal
        gns = co.penatration
        f[co.slave.dofs] -= penalty*gns*N
        #K[co.slave.dofs,co.slave.dofs] += penalty*[co.normal[1], co.normal[2]]*[co.normal[1], co.normal[2]]';
    end

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

