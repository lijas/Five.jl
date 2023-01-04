export LinearSolidElement
"""
LinearSolidElement

"""

struct LinearSolidElement{
        dim,
        order,
        shape,
        T         <: AbstractFloat,
        CV        <: Ferrite.Values,
        DIM       <: MaterialModels.AbstractDim
    } <: AbstractElement{dim}

    thickness::T #used in 2d

    celltype::Type{<:Cell}
    cv::CV
    field::Field
    dimstate::DIM
end

initial_element_state(::LinearSolidElement) = EmptyElementState()

getquadraturerule(e::LinearSolidElement) = e.cv.qr
Ferrite.getnquadpoints(e::LinearSolidElement) = getnquadpoints(e.cv)
Ferrite.getcelltype(e::LinearSolidElement) = e.celltype
Ferrite.ndofs(e::LinearSolidElement) = getnbasefunctions(e.cv)
has_constant_massmatrix(::LinearSolidElement) = true
get_fields(e::LinearSolidElement) = return [e.field]

function LinearSolidElement{dim, order, refshape, T}(;
        thickness = 1.0, 
        qr_order::Int=2, 
        celltype::Type{<:Cell}, 
        dimstate::AbstractDim{dim} = MaterialModels.Dim{3}()
        ) where {dim, order, refshape, T}
    
    ip = Lagrange{dim, refshape, order}()
    geom_ip = Ferrite.default_interpolation(celltype)

    qr = QuadratureRule{dim, refshape}(qr_order)

    cv = CellVectorValues(qr, ip, geom_ip)
    return LinearSolidElement{dim, order, refshape, T, typeof(cv), typeof(dimstate)}(thickness, celltype, cv, Field(:u, ip, dim), dimstate)
end


function integrate_forcevector!(
        element::LinearSolidElement{dim, order, shape, T}, 
        elementstate::AbstractVector{<:AbstractElementState}, 
        material::AbstractMaterial, 
        materialstate::AbstractVector{<:AbstractMaterialState},
        stresses::Vector{<:SymmetricTensor{2,3,T}},
        strains::Vector{<:SymmetricTensor{2,3,T}},
        fe::Vector, 
        cell, 
        Δue::Vector,
        ue::Vector,
        due::Vector,
        dt::T) where {dim, order, shape, T}

    cellvalues = element.cv
  
    reinit!(cellvalues, cell)

    n_basefuncs = getnbasefunctions(cellvalues)

    for q_point in 1:getnquadpoints(cellvalues)

        ∇u = function_gradient(cellvalues, q_point, ue)
        ɛ = symmetric(∇u)

        σ, ∂σ∂ɛ, new_matstate = material_response(element.dimstate, material, ɛ, materialstate[q_point])
        materialstate[q_point] = new_matstate
        stresses[q_point] = dim==3 ? σ : MaterialModels.increase_dim(σ)
        strains[q_point] = dim==3 ? ɛ : MaterialModels.increase_dim(ɛ)
        dΩ = getdetJdV(cellvalues, q_point) * element.thickness
        
        for i in 1:n_basefuncs
            δɛ = shape_symmetric_gradient(cellvalues, q_point, i)
            fe[i] += (σ ⊡ δɛ) * dΩ
        end
    end
    
end

function integrate_forcevector_and_stiffnessmatrix!(
        element::LinearSolidElement{dim, order, shape, T}, 
        elementstate::AbstractVector{<:AbstractElementState}, 
        material::AbstractMaterial, 
        materialstate::AbstractVector{<:AbstractMaterialState},
        stresses::Vector{<:SymmetricTensor{2,3,T}},
        strains::Vector{<:SymmetricTensor{2,3,T}},
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

    for q_point in 1:getnquadpoints(cellvalues)

        ∇u = function_gradient(cellvalues, q_point, ue)
        ɛ = symmetric(∇u)

        σ, ∂σ∂ɛ, new_matstate = material_response(element.dimstate, material, ɛ, materialstate[q_point])
        materialstate[q_point] = new_matstate
        stresses[q_point] = dim==3 ? σ : MaterialModels.increase_dim(σ)
        strains[q_point] = dim==3 ? ɛ : MaterialModels.increase_dim(ɛ)


        dΩ = getdetJdV(cellvalues, q_point) * element.thickness

        for i in 1:n_basefuncs
            δɛi = shape_symmetric_gradient(cellvalues, q_point, i)
            
            fe[i] += (σ ⊡ δɛi) * dΩ
            
            ɛC = δɛi ⊡ ∂σ∂ɛ
            for j in 1:n_basefuncs 
                δɛj = shape_symmetric_gradient(cellvalues, q_point, j)
                ke[i, j] += (ɛC ⊡ δɛj) * dΩ 
            end
        end
    end
    
end

function integrate_massmatrix!(element::LinearSolidElement{dim, order, shape, T, M}, ::AbstractVector{<:AbstractElementState}, material::AbstractMaterial, cell, me::Matrix, ue::AbstractVector, due::AbstractVector) where {dim, order, shape, T, M}

    cv = element.cv
    reinit!(cv, cell)
    ndofs = Ferrite.ndofs(element)

    for qp in 1:getnquadpoints(cv)
        dV = getdetJdV(cv, qp) * element.thickness
        for i in 1:ndofs
            Ni = shape_value(cv, qp, i)
            for j in 1:ndofs
                Nj = shape_value(cv, qp, j)
                me[i, j] += density(material)*Ni⋅Nj*dV;
            end
        end
    end
end


function integrate_dissipation!(
    element       :: LinearSolidElement{dim_p,dim_s,CV},
    elementstate  :: Vector{<:AbstractElementState}, 
    material      :: AbstractMaterial,
    materialstate :: AbstractArray{<:AbstractMaterialState},
    fe            :: Vector{T},
    ge            :: Base.RefValue{T},
    coords, 
    Δue           :: Vector,
    ue            :: Vector,
    due           :: Vector,
    Δt            :: T
    ) where {dim_p,dim_s,CV,T}

    error("function not implemented for $(typeof(element))")
    
end