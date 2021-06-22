"""
Rigid element

"""

struct SpringInfo{dim, T}
    l0::T
    d0::Vec{dim,T}
    #x0_1::Vec{dim,T}
    #x0_2::Vec{dim,T}
    dofs1::Vec{dim,Int}
    dofs2::Vec{dim,Int}
    function SpringInfo{dim,T}(l0::T, d0::Vec{dim,T}) where {dim,T}
        _123 = Vec{dim,T}(Tuple(1:dim))
        _456 = Vec{dim,T}(Tuple(_123.+dim))
        return new(l0, d0, _123, _456)
    end
end

struct SpringElement{dim} <: AbstractElement
    info::Info

    fields::Vector{Field}
    celltype::Type{<:Cell}
    #bcvalue::BCValues{Float64}
end

function SpringElement{dim}() where dim
    return SpringElement{dim}([Field(:u, Lagrange{1, RefCube, 1}(), dim)], Cell{dim,2,2})
end

#Ferrite.ndofs(::SpringElement{dim}) where {dim} = dim*2

function integrate_massmatrix!(element::SpringElement, material::AbstractMaterial, cell::CellIterator, me::Matrix)#, ue::Vector)
    #The spring should always be connected to other parts of the structure, so it does not contribute with mass
end

function integrate_forcevector!(element::SpringElement{dim}, material::AbstractMaterial, fe::Vector, cell::CellIterator, ue::Vector, due::Vector) where {dim,T}
    d = info.d0*info.l0 - (ue[info.dofs1] + ue[info.dofs2])
    
    l = norm(d)
    if norm(d) == 0
        d = Vec{dim,T}((0.0,1.0))
    else
        d /= l
    end

    F = (info.l0 - l)*material.k
    fe[info.dofs1] += F*d
    fe[info.dofs2] -= F*d
end

function integrate_forcevector_and_stiffnessmatrix!(element::SpringElement, material::AbstractMaterial, fe::Vector, ke::AbstractMatrix, cell, ue::Vector, due::Vector)

end

function bodyforce!(element::SpringElement{dim}, material::AbstractMaterial, cell::CellIterator, fe::Vector, forcevec::Vector) where dim

end
