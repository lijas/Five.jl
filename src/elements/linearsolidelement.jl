export LinearSolidElement
"""
LinearSolidElement

"""

struct LinearSolidElement{sdim,CV<:Ferrite.AbstractCellValues, DIMSTATE<:MaterialModels.AbstractDim} <: AbstractElement{sdim}
    celltype::Type{<:Ferrite.AbstractCell}
    cv::CV
    dimstate::DIMSTATE
    thickness::Float64
end

getquadraturerule(e::LinearSolidElement) = getquadraturerule(e.cv)
Ferrite.getnquadpoints(e::LinearSolidElement) = getnquadpoints(e.cv)
Ferrite.getcelltype(e::LinearSolidElement) = e.celltype
Ferrite.ndofs(e::LinearSolidElement) = getnbasefunctions(e.cv)
has_constant_massmatrix(::LinearSolidElement) = true
get_fields(e::LinearSolidElement) = return [(:u, e.cv.ip)]

function LinearSolidElement(;
        celltype::Type{<:Ferrite.AbstractCell{refshape}}, 
        thickness                   = 1.0, 
        qr_order::Int               = 2, 
        dimstate::AbstractDim{sdim} = MaterialModels.Dim{3}(),
        ip::Interpolation           = Ferrite.default_interpolation(celltype),
    ) where {refshape, sdim}
    
    @assert Ferrite.getrefshape(ip) == refshape "refshape of element does not match that of the given interpolation"

    geo_ip = Ferrite.default_interpolation(celltype)
    qr = QuadratureRule{refshape}(qr_order)
    cv = CellValues(qr, ip^sdim, geo_ip^sdim)

    return LinearSolidElement{sdim}(celltype, cv, dimstate, thickness)
end

function integrate_forcevector!(
        element::LinearSolidElement{dim}, 
        elementstate::AbstractVector{<:AbstractElementState}, 
        material::AbstractMaterial, 
        materialstate::AbstractVector{<:AbstractMaterialState},
        stresses::Vector{<:SymmetricTensor{2,3}},
        strains::Vector{<:SymmetricTensor{2,3}},
        fe::Vector, 
        coords, 
        Δue::Vector,
        ue::Vector,
        due::Vector,
        dt::AbstractFloat) where dim

    cellvalues = element.cv
  
    reinit!(cellvalues, coords)

    n_basefuncs = getnbasefunctions(cellvalues)

    for q_point in 1:getnquadpoints(cellvalues)

        ∇u = function_gradient(cellvalues, q_point, ue)
        ɛ = symmetric(∇u)

        σ, ∂σ∂ɛ, new_matstate = material_response(element.dimstate, material, ɛ, materialstate[q_point])
        materialstate[q_point] = new_matstate
        stresses[q_point] = dim==3 ? σ : MaterialModels.increase_dim(σ)
        strains[q_point]  = dim==3 ? ɛ : MaterialModels.increase_dim(ɛ)
        dΩ = getdetJdV(cellvalues, q_point) * element.thickness
        
        for i in 1:n_basefuncs
            δɛ = shape_symmetric_gradient(cellvalues, q_point, i)
            fe[i] += (σ ⊡ δɛ) * dΩ
        end
    end
    
end

function integrate_forcevector_and_stiffnessmatrix!(
        element::LinearSolidElement{dim}, 
        elementstate::AbstractVector{<:AbstractElementState}, 
        material::AbstractMaterial, 
        materialstate::AbstractVector{<:AbstractMaterialState},
        stresses::Vector{<:SymmetricTensor{2,3}},
        strains::Vector{<:SymmetricTensor{2,3}},
        ke::AbstractMatrix, 
        fe::Vector, 
        coords, 
        Δue::Vector,
        ue::Vector,
        due::Vector,
        dt::AbstractFloat) where {dim}

    cellvalues = element.cv
  
    reinit!(cellvalues, coords)

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
            δɛi = symmetric(shape_gradient(cellvalues, q_point, i))
            
            fe[i] += (σ ⊡ δɛi) * dΩ
            
            ɛC = δɛi ⊡ ∂σ∂ɛ
            for j in 1:n_basefuncs 
                δɛj = symmetric(shape_gradient(cellvalues, q_point, j))
                ke[i, j] += (ɛC ⊡ δɛj) * dΩ 
            end
        end
    end
    
end

function integrate_massmatrix!(element::LinearSolidElement{dim}, ::AbstractVector{<:AbstractElementState}, material::AbstractMaterial, cell, me::Matrix, ue::AbstractVector, due::AbstractVector) where {dim}

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
    element       :: LinearSolidElement,
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
    ) where T

    error("function not implemented for $(typeof(element))")
    
end