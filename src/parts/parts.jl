
Ferrite.getdim(p::P) where P<:AbstractPart{dim} where dim = dim

"""
    get_field(::AbstractPart)

Should return the fields of the part.

Example: 
    Solid-element: Field(:u, Lagrange{2,RefCube,1}, 2)
    Shell-element: Field(:u, Lagrange{2,RefCube,1}, 3), Field(:θ, Lagrange{2,RefCube,1}, 3)

"""
get_fields

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
bigtest!(dh::AbstractDofHandler, part::AbstractPart, state::StateVariables)
"""
function bigtest()

end

include("fepart.jl")
include("cohesive_part.jl")
#include("rigidbody/rigidpart.jl")

#init!
#assemble_stiffnessmatrix_and_forcevector!
#get_vtk_part_grid
#get_vtk_part_point_data
