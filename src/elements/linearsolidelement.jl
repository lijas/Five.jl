
"""
LinearSolidElement

Any order, any shape, any dim, 

"""

struct LinearSolidElement{dim,order,shape,T,CV<:Ferrite.Values} <: AbstractElement
    thickness::T #used in 2d

    celltype::Type{<:Cell}
    cv::CV
    field::Field
end


Ferrite.getnquadpoints(e::LinearSolidElement) = getnquadpoints(e.cv)
Ferrite.ndofs(e::LinearSolidElement) = getnbasefunctions(e.cv)
Ferrite.getcelltype(e::LinearSolidElement) = e.celltype
getncoords(s::LinearSolidElement) = Ferrite.getngeobasefunctions(s.cv)

has_constant_massmatrix(::LinearSolidElement) = true

get_fields(e::LinearSolidElement) = return [e.field]

function LinearSolidElement{dim, order, refshape, T}(; thickness = 1.0, qr_order::Int=2, celltype::Type{<:Cell}) where {dim, order, refshape, T}
    
    ip = Lagrange{dim, refshape, order}()
    geom_ip = Ferrite.default_interpolation(celltype)

    qr = QuadratureRule{dim, refshape}(qr_order)

    cv = CellVectorValues(qr, ip, geom_ip)
    return LinearSolidElement{dim, order, refshape, T, typeof(cv)}(thickness, celltype, cv, Field(:u, ip, dim))
end

function integrate_forcevector_and_stiffnessmatrix!(element::LinearSolidElement{dim, order, shape, T}, 
    elementstate::AbstractElementState, 
    material::AbstractMaterial, 
    materialstate::AbstractVector{<:AbstractMaterialState},
    ke::AbstractMatrix, 
    fe::Vector, 
    cell, 
    Δue::Vector,
    ue::Vector,
    due::Vector,
    dt::T) where {dim, order, shape, T}

    cellvalues = element.cv
  
    reinit!(cellvalues, cell)

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
            for j in 1:n_basefuncs 
                ke[i, j] += (ɛC ⊡ δɛ[j]) * dΩ 
            end
        end
    end
    
    return V
end
