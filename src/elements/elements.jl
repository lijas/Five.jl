"""

"""
abstract type AbstractElement end
abstract type AbstractElementState end   

#Many elementstates will be empty, so create default empty one
struct EmptyElementState <: AbstractElementState end

EmptyElementState(::AbstractElement) = EmptyElementState() 
elementstate_type(::Type{T}) where {T<:AbstractElement} = EmptyElementState

"""
    getnquadpoints(::AbstractElement)

Returns the number of quadrature points in the element
Used for allocating arrays for the material states
"""
getnquadpoints

"""
    getquadraturerule(::AbstractElement)

Returns the number of quadrature rule used by the element
Used for postprocessing and L2Projection
"""
getquadraturerule

"""
    get_fields(::AbstractElement)

Returns a list with fields of element, e.g :u, :T, :θ 
"""
function get_fields(::AbstractElement) end

"""
    ndofs(::AbstractElement)

Returns the number of dofs in the element
Used for pre-allocating element matrices, internal force matrices etc.
"""
ndofs

"""
    getcelltype(::AbstractElement)

Returns the celltype (Ferrite.AbstractCell) of the element
"""
getcelltype


"""
    Returns true if the massmatrix is constant or not
"""
has_constant_massmatrix

"""
    Returns the internal forcevector of element
"""
integrate_forcevector!

"""
    Returns the mass matrix of the element
"""
integrate_massmatrix!
 
"""
    Returns the internal forcevector and the element stiffness matrix of element
"""
integrate_forcevector_and_stiffnessmatrix!

"""
    Returns the dissipation and the gradient of the dissipation (used in LocalDissipationSolver)
"""
integrate_dissipation!

"""
    Returns the fstar vector, used in the DissipationSolver
"""
integrate_fstar!

include("solidelement.jl")
include("bar_element.jl")
include("linearsolidelement.jl")

#Cohesive
include("../utils/cohesive_element_utils.jl")
include("cohesive_element.jl")



# # # # #
# Element defintions
# # # # #
export SolidElementQuad, SolidElementHexa, SolidElementTria, SolidElementTetra

const LinearSolidElementQuad = LinearSolidElement{2,1,RefCube,Float64}
const LinearSolidElementTria = LinearSolidElement{2,1,RefTetrahedron,Float64}
const LinearSolidElementHexa = LinearSolidElement{3,1,RefCube,Float64}

const SolidElementQuad = SolidElement{2,1,RefCube,Float64}
const SolidElementHexa = SolidElement{3,1,RefCube,Float64}
const SolidElementTria = SolidElement{2,1,RefTetrahedron,Float64}
const SolidElementTetra = SolidElement{3,1,RefTetrahedron,Float64}