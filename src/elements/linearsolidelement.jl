
"""
LinearSolidElement

Any order, any shape, any dim, 

"""

struct LinearSolidElement{dim,order,shape,T,CV<:JuAFEM.Values} <: AbstractElement

    thickness::T #used in 2d

    ndofs::Int
    fields::Vector{Field}
    celltype::Type{<:Cell}
    M2::Int
    cv::CV
end

struct LinearSolidElementState <: AbstractElementState
end

LinearSolidElementState(::LinearSolidElement{dim,order,shape,T}, 
                        ::AbstractVector{Vec{dim,T}} = Vec{dim,T}[]) where {dim,order,shape,T} = LinearSolidElementState()

get_elementstate_type(e::LinearSolidElement) = LinearSolidElementState
JuAFEM.getnquadpoints(e::LinearSolidElement) = JuAFEM.getnquadpoints(e.cv)
JuAFEM.ndofs(e::LinearSolidElement) = e.ndofs
JuAFEM.getcelltype(e::LinearSolidElement) = e.celltype
getfaceinterpolation(element::LinearSolidElement{2,1,shape}) where {shape} = Lagrange{1, RefCube, 1}()
getfaceinterpolation(element::LinearSolidElement{3,1,shape}) where {shape} = Lagrange{2, shape, 1}()
has_constant_massmatrix(::LinearSolidElement) = true
getncoords(s::LinearSolidElement) = JuAFEM.getngeobasefunctions(s.cv)

function LinearSolidElement{dim, order, refshape, T}(; thickness::T = 1.0) where {dim, order, refshape, T}
    
    ip = Lagrange{dim, refshape, order}()
    qr = QuadratureRule{dim, refshape}(2)
    nnodes = getnbasefunctions(ip)
    nfaces = length(JuAFEM.faces(ip))
    M2 =  dim == 2 ? 3 : 6
    ndofs = nnodes*dim
    cv = CellVectorValues(qr, ip)
    return LinearSolidElement{dim, order, refshape, T, typeof(cv)}(thickness, ndofs, [Field(:u, ip, dim)], Cell{dim,nnodes,nfaces}, M2, cv)
end

function LinearSolidElement{dim, order, refshape, T}(cv::JuAFEM.Values{dim}, ip) where {dim, order, refshape, T}
    
    nnodes = JuAFEM.getngeobasefunctions(cv)
    nfaces = length(JuAFEM.faces(ip))
    M2 =  dim == 2 ? 3 : 6
    ndofs = JuAFEM.getnbasefunctions(cv)
    return LinearSolidElement{dim, order, refshape, T, typeof(cv)}(ndofs, [Field(:u, ip, dim)], Cell{dim,nnodes,nfaces}, M2, cv)
end

function integrate_forcevector_and_stiffnessmatrix!(element::LinearSolidElement{dim, order, shape, T, M}, 
    elementstate::LinearSolidElementState, 
    material::AbstractMaterial, 
    materialstate::AbstractVector{<:AbstractMaterialState},
    ke::AbstractMatrix, 
    fe::Vector, 
    cell, 
    Δue::Vector,
    ue::Vector,
    due::Vector,
    dt::T) where {dim, order, shape, T, M}

    cellvalues = element.cv
  
    @timeit "reinit" reinit!(cellvalues, cell)

    n_basefuncs = getnbasefunctions(cellvalues)

    δɛ = [zero(SymmetricTensor{2, dim}) for i in 1:n_basefuncs]
    V = 0
    for q_point in 1:getnquadpoints(cellvalues)

        ∇u = function_gradient(cellvalues, q_point, ue)
        ɛ = symmetric(1/2 * (∇u + ∇u'))
   
        #MatLinearElasticState(zero(SymmetricTensor{2,dim,T}))
        σ, ∂σ∂ɛ, new_matstate = constitutive_driver(material, ɛ, materialstate[q_point])
        materialstate[q_point] = new_matstate

        for i in 1:n_basefuncs
            δɛ[i] = symmetric(shape_gradient(cellvalues, q_point, i)) 
        end
        dΩ = getdetJdV(cellvalues, q_point) * element.thickness
        V += dΩ

        for i in 1:n_basefuncs
            δu = shape_value(cellvalues, q_point, i)
            
            fe[i] += (σ ⊡ δɛ[i]) * dΩ
            
            ɛC = δɛ[i] ⊡ ∂σ∂ɛ
            for j in 1:n_basefuncs # assemble only upper half
                ke[i, j] += (ɛC ⊡ δɛ[j]) * dΩ # can only assign to parent of the Symmetric wrapper
            end
        end
    end
    
    return V
end
