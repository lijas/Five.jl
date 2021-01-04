struct CellIterator2{dim,T}
    flags::UpdateFlags
    current_cellid::JuAFEM.ScalarWrapper{Int}
    nodes::Vector{Int}
    coords::Vector{Vec{dim,T}}
    dh::JuAFEM.AbstractDofHandler
    celldofs::Vector{Int}
    cellset::Vector{Int}

    function CellIterator2{dim,T}(dh::JuAFEM.AbstractDofHandler, cellset::Vector{Int}, flags::UpdateFlags) where {dim,T}
        N = length(dh.grid.cells[first(cellset)].nodes)#nnodes_per_cell(dh, first(cellset))
        cell = JuAFEM.ScalarWrapper(0)
        nodes = zeros(Int, N)
        coords = zeros(Vec{dim,T}, N)
        n = ndofs_per_cell(dh, first(cellset))
        celldofs = zeros(Int, n)
        return new{dim,T}(flags, cell, nodes, coords, dh, celldofs, cellset)
    end

end

function CellIterator2(dh::JuAFEM.AbstractDofHandler, element::AbstractElement, cellset::Vector{Int}, flags::UpdateFlags=UpdateFlags()) 
    dim = JuAFEM.getdim(element)
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
@inline JuAFEM.getnodes(ci::CellIterator2) = ci.nodes
@inline JuAFEM.getcoordinates(ci::CellIterator2) = ci.coords
@inline JuAFEM.onboundary(ci::CellIterator2, face::Int) = ci.grid.boundary_matrix[face, ci.current_cellid[]]
@inline JuAFEM.cellid(ci::CellIterator2) = ci.current_cellid[]
@inline JuAFEM.celldofs!(v::Vector, ci::CellIterator2) = celldofs!(v, ci.dh, ci.current_cellid[])
@inline JuAFEM.celldofs(ci::CellIterator2) = ci.celldofs

function JuAFEM.reinit!(ci::CellIterator2{dim}, i::Int) where {dim}
    ci.current_cellid[] = ci.cellset[i]
    
    ci.flags.nodes  && JuAFEM.cellnodes!(ci.nodes, ci.dh, ci.current_cellid[])
    ci.flags.coords && JuAFEM.cellcoords!(ci.coords, ci.dh, ci.current_cellid[])
    ci.flags.celldofs && JuAFEM.celldofs!(ci.celldofs, ci)

    return ci
end
