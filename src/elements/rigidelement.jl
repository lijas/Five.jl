"""
Rigid element

"""
struct RigidElement{dim} <: AbstractElement
    ndofs::Int
    fields::Vector{Field}
    celltype::Type{<:Cell}
    #x::RigidCoordinate
end

function RigidElement{dim}() where dim
    ndofs = dim == 2 ? 3 : 7
    fields = [Field(:xyθ, Serendipity{dim,RefCube,0}(), ndofs)]#, 
              #Field(:λ, Serendipity{dim,RefCube,0}(), 1)]
    return RigidElement{dim}(ndofs, fields, Cell{dim,1,0})
end

struct RigidInfo{dim, T}
    mass::T
    inertia::Tensor{2, 3, T, 9}       #In 2d, the first element is used 
    initial_position::Vec{dim, T} 
    initial_rotation::NTuple{4,T} #In 2d, the first element is used 
    nodes::Vector{Int} #Should probobly be coordinates, node ints...
end

get_elementinfo_type(e::RigidElement{dim}) where {dim} = RigidElement{dim,Float64}
has_constant_massmatrix(::RigidElement{2}) = true
has_constant_massmatrix(::RigidElement{3}) = false

function integrate_massmatrix!(element::RigidElement{2}, state::RigidElementState, info::RigidInfo, material::AbstractMaterial, cell, me::Matrix, ue::AbstractVector, due::AbstractVector)
    me[1,1] = info.mass
    me[2,2] = info.mass
    me[3,3] = info.inertia[1,1]
end

const rigid_alpha = 100
function integrate_massmatrix!(element::RigidElement{3}, state::RigidElementState, info::RigidInfo, material::AbstractMaterial, cell::CellIterator, me::Matrix, ue::AbstractVector, due::AbstractVector)
    me[1,1] = info.mass
    me[2,2] = info.mass
    me[3,3] = info.mass

    θ = info.initial_rotation .+ ue[4:7]
    G = getG(θ)

    C = θ'*θ-1
    Cq = 2*θ'
    me[[4,5,6,7],[4,5,6,7]] = G'*info.inertia*G + 10000*Cq'*Cq

    return me
end

#=
bodydofs_rot = getbodydofs_rot(dh, bodyid)
        bodydofs_trans = getbodydofs_trans(dh, bodyid)

        G = getG(qq[bodyid].θ)
        mtt = G'*body.I*G

        M[bodydofs_trans, bodydofs_trans] = body.M
        M[bodydofs_rot, bodydofs_rot] = mtt

        wbar = G*dqq[bodyid].θ

        #Remember, in 2d dG is equal to zero so this is wrong...
        dG = getG(dqq[bodyid].θ)

        Qv[bodydofs_rot] = -2*wbar'*body.I*dG;
=#

function integrate_forcevector!(element::RigidElement{2}, state::RigidElementState, elementinfo::RigidInfo, material::AbstractMaterial, mstate::Vector{<:AbstractMaterialState}, fe::Vector, cell, ue::Vector, due::Vector)

end

function integrate_forcevector!(element::RigidElement{3}, state::RigidElementState, info::RigidInfo, material::AbstractMaterial, mstate::Vector{<:AbstractMaterialState}, fe::Vector, cell, ue::Vector, due::Vector)
    #Can the quadratic veloctiy vector be seen as internal forces??
    θ = info.initial_rotation .+ ue[4:7]
    dθ = due[4:7]

    G = getGbar(θ)
    wbar = G*dθ 
    dG = getGbar(dθ)

    C = θ'*θ-1
    dC = dθ'*dθ-1

    Cq = 2*θ'

    dCq = 2*dθ'

    fe[4:7] = (-2*wbar'*info.inertia*dG)' +10000*Cq'*(2*dθ'*dθ)# + 2*rigid_alpha*rigid_alpha*dC + rigid_alpha^2*C)
end

function integrate_forcevector_and_stiffnessmatrix!(element::RigidElement{2}, state::RigidElementState, elementinfo::RigidInfo, material::AbstractMaterial, mstate::Vector{AbstractMaterialState}, fe::Vector, ke::AbstractMatrix, cell, ue::Vector, due::Vector)

end

function bodyforce!(element::RigidElement{dim}, state::RigidElementState, info::RigidInfo, material::AbstractMaterial, cell, fe::Vector, forcevec::Vector) where dim
    for i in 1:dim
        fe[i] = info.mass*forcevec[i]
    end
end

function export_vtk_displacements!(data, dh::DofHandler{dim,T}, element::RigidElement{dim}, info::RigidInfo{dim,T}, cell, ue::AbstractVector{T}) where {dim,T}
    
    #ToDo: some way of exporting rigid bodies to vtk
    for i in 1:dim
        data[i,cell.nodes[1]] = ue[i]
    end
    return

end