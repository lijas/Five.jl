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


function Ferrite.close!(ef::Constraints, dh::MixedDofHandler)
    for (i, e) in enumerate(ef.external_forces)
        ForceType = typeof(e)
        ef.external_forces[i] = init_external_force!(e, dh)
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
        apply_external_force!(ef, state, globaldata)
    end
end

#=function _applypenalty!(consteq::DofConstraint, f, K, penalty::T, ue) where T
    @timeit "f" g::T = consteq.func(ue)
    @timeit "G" G::Array{T,1} = consteq.grad(ue)
    @timeit "H" H::Array{T,2} = consteq.hess(ue)*g
    
    f[consteq.dofs] -= penalty*g*G;
    #K[consteq.dofs, consteq.dofs] -= (penalty*G*G' + penalty*H)
end=#

export TestConstraint

struct TestConstraint{BI <: Ferrite.BoundaryIndex}

    faces::Set{BI}
    mastervertex::VertexIndex
    dofs::Vector{Int}
    field::Symbol
end

function TestConstraint(; faces, mastervertex, field, dofs)
    TestConstraint(faces, mastervertex, dofs, field)
end

function create_linear_constraints(dh::MixedDofHandler, tc::TestConstraint{BI}) where BI

    fh = getfieldhandler(dh, cellid(first(tc.faces)))

    masterdofs = dofs_on_vertex(dh, fh, tc.mastervertex, tc.field, tc.dofs);
    @assert(length(masterdofs) == 1)
    masterdof = masterdofs[1]
    slavedofs = Int[]
    for face in tc.faces
        dofs = []
        @show typeof(face)
        if face isa FaceIndex
            dofs = dofs_on_face(dh, fh, face, tc.field, tc.dofs)
        elseif face isa VertexIndex
            dofs = dofs_on_vertex(dh, fh, face, tc.field, tc.dofs)
        else
            error("Hej")
        end
        append!(slavedofs, dofs)
    end

    #@assert( masterdof in slavedofs)
    unique!(slavedofs)
    filter!((x)->x!=masterdof, slavedofs)


    lcs = Ferrite.LinearConstraint[]
    for d in slavedofs
        @show Ferrite.LinearConstraint(d, [masterdof=>1.0], 0.0)
        push!(lcs, Ferrite.LinearConstraint(d, [masterdof=>1.0], 0.0))
    end
     
    return lcs;
end




