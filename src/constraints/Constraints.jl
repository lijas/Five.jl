"""
Constraints
    WIP
"""

abstract type AbstractConstraint end
abstract type AbstractConstraintInfo end

#=struct Constraints{DH} #ContraintHandler
    penalty_constraints::Vector{AbstractConstraint}
    #lagrange_constraints::Vector{AbstractConstraint}

    constraint_infos::Vector{AbstractConstraintInfo}

    penalty_factors::Vector{Float64}

    constraint_dofs::Vector{Int}
    constraint_dofs_offset::Vector{Int}

    dh::DH
end=#

struct Constraints
    external_forces::Vector{AbstractExternalForce}
end

@inline nconstraints(c::Constraints) = length(c.penalty_constraints) + length(c.lagrange_constraints)

#=function Constraints(dh::DH) where {DH}
    Constraints(AbstractConstraint[], AbstractConstraintInfo[], Float64[], Int[], Int[], dh)
end=#

function Constraints()
    Constraints(AbstractExternalForce[])
end

function add_penalty_constraints!(c::Constraints, constraint::AbstractConstraint, penalty_factor::T) where T
    push!(c.penalty_constraints, constraint)
    push!(c.penalty_factors, penalty_factor)
end

function add_lagrange_constraints!(c::Constraints, constraint::AbstractConstraint)
    push!(c.lagrange_constraints, constraint)
end

function close_ch!(ch::Constraints)

    push!(ch.constraint_dofs_offset,1)
    for pc in ch.penalty_constraints

        info, dofs = construct_contraint(pc,ch.dh)
        
        append!(ch.constraint_dofs, dofs)
        push!(ch.constraint_dofs_offset, length(ch.constraint_dofs)+1)
        push!(ch.constraint_infos, info)
    end

end

function constraintdofs(ch, i::Int)
    return ch.constraint_dofs[(ch.constraint_dofs_offset[i]):(ch.constraint_dofs_offset[i+1]-1)]
end

#=function apply_constraints!(ch, f, M, u)
    for (i, c) in enumerate(ch.penalty_constraints)
        dofs = constraintdofs(ch,i)
        ue = u[dofs]
        
        g,G = get_g_and_G(c, ch.constraint_infos[i], ue)

        f[dofs] -= ch.penalty_factors[i]*G*g;
        #M[dofs,dofs] += ch.penalty_factors[i]*G*G';
    end
    #for (i, constraint) in enumerate(c.lagrange_constraints)
        #ue = u[constraint.dofs]
        #_applylagrange!(constraint, f c.penalty_factors[i], ue)
    #end
end=#

function apply_constraints!(dh, ch::Constraints, state::StateVariables, globaldata)
    for ef in ch.external_forces
        _apply_external_force!(dh, ef, state, globaldata)
    end
end

#=function _applypenalty!(consteq::DofConstraint, f, K, penalty::T, ue) where T
    @timeit "f" g::T = consteq.func(ue)
    @timeit "G" G::Array{T,1} = consteq.grad(ue)
    @timeit "H" H::Array{T,2} = consteq.hess(ue)*g
    
    f[consteq.dofs] -= penalty*g*G;
    #K[consteq.dofs, consteq.dofs] -= (penalty*G*G' + penalty*H)
end=#