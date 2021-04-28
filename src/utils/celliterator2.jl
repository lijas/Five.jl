struct CellIterator2{dim,T}
    flags::UpdateFlags
    current_cellid::Ferrite.ScalarWrapper{Int}
    nodes::Vector{Int}
    coords::Vector{Vec{dim,T}}
    dh::Ferrite.AbstractDofHandler
    celldofs::Vector{Int}
    cellset::Vector{Int}

    function CellIterator2{dim,T}(dh::Ferrite.AbstractDofHandler, cellset::Vector{Int}, flags::UpdateFlags) where {dim,T}
        N = length(dh.grid.cells[first(cellset)].nodes)#nnodes_per_cell(dh, first(cellset))
        cell = Ferrite.ScalarWrapper(0)
        nodes = zeros(Int, N)
        coords = zeros(Vec{dim,T}, N)
        n = ndofs_per_cell(dh, first(cellset))
        celldofs = zeros(Int, n)
        return new{dim,T}(flags, cell, nodes, coords, dh, celldofs, cellset)
    end

end

function CellIterator2(dh::Ferrite.AbstractDofHandler, element::AbstractElement, cellset::Vector{Int}, flags::UpdateFlags=UpdateFlags()) 
    dim = Ferrite.getdim(element)
    T = Float64
    return CellIterator2{dim,T}(dh, cellset, flags)
end

# iterator interface
function Base.iterate(ci::CellIterator2, state = 1)
    if state > length(ci.cellset)
        return nothing
    else
        return (reinit!(ci, state), state+1)
    end
end
Base.length(ci::CellIterator2)  = length(ci.cellset)

Base.IteratorSize(::Type{T})   where {T<:CellIterator2} = Base.HasLength() # this is default in Base
Base.IteratorEltype(::Type{T}) where {T<:CellIterator2} = Base.HasEltype() # this is default in Base
Base.eltype(::Type{T})         where {T<:CellIterator2} = T

# utility
@inline Ferrite.getnodes(ci::CellIterator2) = ci.nodes
@inline Ferrite.getcoordinates(ci::CellIterator2) = ci.coords
@inline Ferrite.onboundary(ci::CellIterator2, face::Int) = ci.grid.boundary_matrix[face, ci.current_cellid[]]
@inline Ferrite.cellid(ci::CellIterator2) = ci.current_cellid[]
@inline Ferrite.celldofs!(v::Vector, ci::CellIterator2) = celldofs!(v, ci.dh, ci.current_cellid[])
@inline Ferrite.celldofs(ci::CellIterator2) = ci.celldofs

function Ferrite.reinit!(ci::CellIterator2{dim}, i::Int) where {dim}
    ci.current_cellid[] = ci.cellset[i]
    
    ci.flags.nodes  && Ferrite.cellnodes!(ci.nodes, ci.dh, ci.current_cellid[])
    ci.flags.coords && Ferrite.cellcoords!(ci.coords, ci.dh, ci.current_cellid[])
    ci.flags.celldofs && Ferrite.celldofs!(ci.celldofs, ci)

    return ci
end
