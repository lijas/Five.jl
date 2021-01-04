
"""
ThroughThicknessStresses
    Struct containing the stresses, strains, and coordinates in the local and global systems
"""
struct ThroughThicknessStresses
    stresses::Vector{SymmetricTensor{2,3,Float64,6}}
    strains::Vector{SymmetricTensor{2,3,Float64,6}}
    local_coords::Vector{Vec{3,Float64}}
    global_coords::Vector{Vec{3,Float64}}
    function ThroughThicknessStresses(stresses::Vector{SymmetricTensor{2,3,T,M}}, 
                                      strains::Vector{SymmetricTensor{2,3,T,M}},
                                      zcoords::Vector{Vec{3,T}},
                                      coords::Vector{Vec{3,T}}) where {T,M}
        @assert length(coords) == length(stresses)
        @assert length(coords) == length(strains)
        @assert length(coords) == length(zcoords)
        new(stresses, strains, zcoords, coords)
    end
end

function ThroughThicknessStresses() 
    ThroughThicknessStresses(SymmetricTensor{2,3,Float64,6}[], SymmetricTensor{2,3,Float64,6}[], Vec{3,Float64}[], Vec{3,T}[])
end