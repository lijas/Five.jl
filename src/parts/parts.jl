
Ferrite.getdim(p::P) where P<:AbstractPart{dim} where dim = dim

"""
    get_field(::AbstractPart)

Should return the fields of the part.

Example: 
    Solid-element: Field(:u, Lagrange{2,RefCube,1}, 2)
    Shell-element: Field(:u, Lagrange{2,RefCube,1}, 3), Field(:θ, Lagrange{2,RefCube,1}, 3)

"""
function get_fields(::AbstractPart) end

"""
    init_part!
"""
init_part!

""" 
    assemble_stiffnessmatrix_and_forcevector!(dh::AbstractDofHandler, part::AbstractPart, state::StateVariables)
"""
assemble_stiffnessmatrix_and_forcevector!

"""
    assemble_forcevector!(dh::AbstractDofHandler, part::AbstractPart, state::StateVariables)
"""
assemble_forcevector!

"""
    assemble_dissipation!(dh::AbstractDofHandler, part::AbstractPart, state::StateVariables)
"""
assemble_dissipation!

"""
    assemble_sparsity_pattern!(dh::AbstractDofHandler, part::AbstractPart, state::StateVariables)
"""
function assemble_sparsity_pattern!(::AbstractPart, globaldata)
    return sparse(Int[], Int[], Float64[], ndofs(globaldata.dh), ndofs(globaldata.dh))
end

