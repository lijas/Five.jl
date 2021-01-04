
JuAFEM.getdim(p::P) where P<:AbstractPart{dim} where dim = dim

include("fepart.jl")
include("igapart.jl")
#include("rigidbody/rigidpart.jl")

#init!
#assemble_stiffnessmatrix_and_forcevector!
#get_vtk_part_grid
#get_vtk_part_point_data
