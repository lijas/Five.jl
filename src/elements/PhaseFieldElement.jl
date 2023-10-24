export PhaseFieldElement
"""
PhaseFieldElement

"""

struct PhaseFieldElement{dim, CV1<:Ferrite.AbstractCellValues,CV2<:Ferrite.AbstractCellValues,DIM<:MaterialModels.AbstractDim} <: AbstractElement{dim}
    celltype::Type{<:Ferrite.AbstractCell}
    cv_u::CV1
    cv_d::CV2
    dimstate::DIM
    thickness::Float64
end

is_dissipative(::PhaseFieldElement) = true

getquadraturerule(e::PhaseFieldElement) = e.cv_u.qr
Ferrite.getnquadpoints(e::PhaseFieldElement) = getnquadpoints(e.cv_u)
Ferrite.ndofs(e::PhaseFieldElement) = getnbasefunctions(e.cv_u) + getnbasefunctions(e.cv_d)
Ferrite.getcelltype(e::PhaseFieldElement) = e.celltype
getncoords(s::PhaseFieldElement) = Ferrite.getngeobasefunctions(s.cv_u)

has_constant_massmatrix(::PhaseFieldElement) = true

get_fields(e::PhaseFieldElement{dim,order,shape,T}) where {dim,order,shape,T} = return [(:u, e.cv_u.ip), (:d, e.cv_d.ip)]

function PhaseFieldElement(;
    celltype::Type{<:Ferrite.AbstractCell},
    thickness         = 1.0,     
    qr_order::Int     = 2,
    ip::Interpolation = Ferrite.default_interpolation(celltype), 
    dimstate::Optional{MaterialModels.AbstractDim} = nothing) where {}
    
    geom_ip = Ferrite.default_interpolation(celltype)
    refshape = Ferrite.getrefshape(geom_ip)
    sdim = Ferrite.getdim(ip)

    qr = QuadratureRule{refshape}(qr_order)
    cv_u = CellValues(qr, ip^sdim, geom_ip^sdim)
    cv_d = CellValues(qr, ip, geom_ip^sdim)

    if dimstate === nothing
        if sdim == 3
            dimstate = MaterialModels.Dim{3}()
        else sdim == 2
            dimstate = MaterialModels.PlaneStrain()
        end
    end

    return PhaseFieldElement(celltype, cv_u, cv_d, dimstate, thickness)
end

function integrate_forcevector_and_stiffnessmatrix!(
        element::PhaseFieldElement{dim,}, 
        elementstate::AbstractVector{<:AbstractElementState}, 
        material::AbstractMaterial, 
        materialstate::AbstractVector{<:AbstractMaterialState},
        stresses::Vector{<:SymmetricTensor{2,3,T}},
        strains::Vector{<:SymmetricTensor{2,3,T}},
        ke::Matrix,
        fe::Vector, 
        cell, 
        Δue::Vector,
        a::Vector,
        due::Vector,
        Δt::T) where {dim, T}

        
    cv_u = element.cv_u
    cv_d = element.cv_d
    reinit!(cv_u, cell)
    reinit!(cv_d, cell)
    ndofs_u = getnbasefunctions(cv_u)
    ndofs_d = getnbasefunctions(cv_d)

    urange = (1:ndofs_u) .+ 0
    drange = (1:ndofs_d) .+ ndofs_u

    lc = material.lc
    Gc = material.Gc

    for qp in 1:getnquadpoints(cv_u)
        
        ∇u = function_gradient(cv_u, qp, a, urange)
        ∇d = function_gradient(cv_d, qp, a, drange)
        d = function_value(cv_d, qp, a, drange)
        ε = symmetric(∇u)
        
        dΩ = getdetJdV(cv_u, qp) * element.thickness

        #σ, dσdε, σ⁺, Ψ⁺ = phasefield_response(ε, d, Gc, λ, μ)
        #σ, dσdε, state = material_response(element.dimstate, material, ε, d, materialstate[qp])
        prev_state = materialstate[qp]
        σ, dσdε, new_state = my_material_response(material, ε, d, materialstate[qp])

        σ⁺ = new_state.σ⁺
        H = new_state.Ψ⁺
        ∂H∂ε = σ⁺
        materialstate[qp] = new_state
        if H <= prev_state.Ψ⁺
            H = prev_state.Ψ⁺
            ∂H∂ε = zero(σ⁺)
            materialstate[qp] = prev_state
        end

        if dim == 2
            ∂H∂ε = MaterialModels.reduce_dim(∂H∂ε, element.dimstate)
            σ⁺   = MaterialModels.reduce_dim(σ⁺  , element.dimstate)
        end

        #Store stress and strain (always as 3d)
        stresses[qp] = dim==2 ? MaterialModels.increase_dim(σ) : σ 
        strains[qp]  = dim==2 ? MaterialModels.increase_dim(ε) : ε
        
        for i in 1:ndofs_u
            ∇Ni = shape_symmetric_gradient(cv_u, qp, i)

            fe[i] += (σ ⊡ ∇Ni) * dΩ

            for j in 1:ndofs_u
                ∇Nj = shape_symmetric_gradient(cv_u, qp, j)

                ke[i,j] += (∇Ni ⊡ dσdε ⊡ ∇Nj) *dΩ 
            end

            for j in 1:ndofs_d
                Nj = shape_value(cv_d, qp, j)
                
                ke[i,j+ndofs_u] += (2*(d-1)*Nj*(∇Ni ⊡ σ⁺)) *dΩ 
            end
        end

        for i in 1:ndofs_d
            Ni = shape_value(cv_d, qp, i)
            ∇Ni = shape_gradient(cv_d, qp, i)

            fe[i+ndofs_u] += (((Gc/lc + 2*H)*d - 2*H)*Ni + Gc*lc*∇d⋅∇Ni) * dΩ

            for j in 1:ndofs_u
                ∇Nj = shape_symmetric_gradient(cv_u, qp, j)

                ke[i+ndofs_u,j] += (2*(d-1)*Ni*(∇Nj ⊡ ∂H∂ε)) *dΩ 
            end

            for j in 1:ndofs_d
                Nj = shape_value(cv_d, qp, j)
                ∇Nj = shape_gradient(cv_d, qp, j)

                ke[i+ndofs_u, j+ndofs_u] += ((Gc/lc + 2*H)*Ni*Nj + Gc*lc*∇Ni⋅∇Nj) *dΩ 
            end

        end

    end

end


function integrate_dissipation!(
    element       :: PhaseFieldElement{dim}, 
    elementstate  ::AbstractVector{<:AbstractElementState}, 
    material      :: AbstractMaterial,
    materialstate :: AbstractArray{<:AbstractMaterialState},
    fe            :: Vector{T},
    ge            :: Base.RefValue{T},
    coords, 
    a             :: Vector,
    Δa            :: Vector,
    due           :: Vector,
    Δt            :: T
    ) where {dim, T}
    
    cv_u = element.cv_u
    cv_d = element.cv_d
    reinit!(cv_d, coords)
    ndofs_u = getnbasefunctions(cv_u)
    ndofs_d = getnbasefunctions(cv_d)

    de = a[(1:ndofs_d) .+ ndofs_u]
    Δde = Δa[(1:ndofs_d) .+ ndofs_u]

    lc = material.lc
    Gc = material.Gc

    for qp in 1:getnquadpoints(cv_u)
        
        ∇d = function_gradient(cv_d, qp, de)
        d = function_value(cv_d, qp, de)

        Δd = function_value(cv_d, qp, Δde)
        Δ∇d = function_gradient(cv_d, qp, Δde)
        
        dΩ = getdetJdV(cv_u, qp) * element.thickness
       
        ge[] += (Gc/lc) * (d⋅Δd + lc^2*∇d⋅Δ∇d) * dΩ
        for i in 1:ndofs_d
            Ni = shape_value(cv_d, qp, i)
            ∇Ni = shape_gradient(cv_d, qp, i)

            fe[i+ndofs_u] += (Gc/lc) * (Ni*Δd + d*Ni + lc^2*(∇Ni⋅Δ∇d + ∇d⋅∇Ni)) * dΩ
        end

    end
end

#mac₊(x::T) where T = x < 0.0 ? 0.0 : x #0.5(x + abs(x))
#mac₋(x::T) where T = x > 0.0 ? 0.0 : x #0.5(x - abs(x))
#heaviside(x::T) where T = x > 0.0 ? one(T) : zero(T)
#=
function phasefield_response_AD(ε::SymmetricTensor{2,dim}, φ, Gc, λ, μ) where {dim}
    ε_dual = Tensors._load(ε, nothing)
    _σ, _σ⁺, _ψ⁺ = _phasefield_response(ε_dual, φ, Gc, λ, μ)

    σ = Tensors._extract_value(_σ)
    dσdε = Tensors._extract_gradient(_σ, ε)
    σ⁺ = Tensors._extract_value(_σ⁺)
    Ψ⁺ = Tensors._extract_value(_ψ⁺)

    return σ, dσdε, σ⁺, Ψ⁺
end

function _phasefield_response2(ε::SymmetricTensor{2,dim,T}, φ, Gc, λ, μ) where {dim,T}
    I = one(SymmetricTensor{2,dim})
    σ₀  = (λ*tr(ε)*I + 2*μ*ε)
    σ  = (1-φ)^2 * σ₀
    ψ⁺ = (0.5*λ*tr(ε)^2 + μ*tr(ε ⋅ ε))
    return σ, σ₀, ψ⁺
end

function _phasefield_response(ε::SymmetricTensor{2,dim,T}, φ, Gc, λ, μ) where {dim,T}
    I = one(SymmetricTensor{2,dim})

    p, n = Tensors.eigen(ε)
    ε⁺ = zero(SymmetricTensor{2,dim,T})
    ε⁻ = zero(SymmetricTensor{2,dim,T})
    
    for i in 1:dim
        ε⁺ += mac₊(p[i]) * symmetric((n[:,i] ⊗ n[:,i]))
        ε⁻ += mac₋(p[i]) * symmetric((n[:,i] ⊗ n[:,i]))
    end

    σ⁺ = λ*mac₊(tr(ε))*I + 2*μ*ε⁺
    σ⁻ = λ*mac₋(tr(ε))*I + 2*μ*ε⁻

    σ = (1-φ)^2 * σ⁺ + σ⁻ 
    ψ⁺ = (0.5*λ*mac₊(tr(ε))^2 + μ*(ε⁺ ⊡ ε⁺))

    return σ, σ⁺, ψ⁺
end


function phasefield_response(ε::SymmetricTensor{2,dim}, φ, Gc, λ, μ) where {dim}

    I2 = one(SymmetricTensor{2,dim})
    I4 = one(SymmetricTensor{4,dim})
    I2oI2 = symmetric(I2⊗I2)

    gφ = (1-φ)^2

    I⁺ = positive_proj(ε)
    I⁻ = I4 - I⁺

    ε⁺ = I⁺ ⊡ ε
    ε⁻ = ε - ε⁺

    trε = tr(ε)

    σ⁺ = λ * mac₊(trε) * I2 + 2.0*Gc * ε⁺
    σ⁻ = λ * mac₋(trε) * I2 + 2.0*Gc * ε⁻
    σ = gφ*σ⁺ + σ⁻

    Ψ⁺ = 0.5*λ*mac₊(trε)^2 + Gc*(ε⁺⊡ε⁺)

    H1 = heaviside(trε)
    H2 = 1.0 - H1

    dσ = gφ * (2.0*Gc*I⁺ + λ*H1*(I2oI2)) + 
              (2.0*Gc*I⁻ + λ*H2*(I2oI2))

    return σ, dσ, σ⁺, Ψ⁺
end

function positive_proj(ε::SymmetricTensor{2,dim}) where dim

    I⁺ = zero(SymmetricTensor{4,dim})
    epos = zeros(dim)
    diag = zeros(dim)

    v, p = eigen(ε)
    for i in 1:dim
        epos[i] = mac₊(v[i])
        if v[i] > 0.0
            diag[i]=1.0
        end
    end

    for i in 1:dim
        Ma = p[:,i] ⊗ p[:,i]
        I⁺ += symmetric(Ma ⊗ Ma) * diag[i]
    end

    for a in 1:dim
        for b in 1:dim
            Ma = p[:,a] ⊗ p[:,a]
            Mb = p[:,b] ⊗ p[:,b]

            Gab = otimesu(Ma,Mb) + otimesl(Ma,Mb)
            Gba = otimesu(Mb,Ma) + otimesl(Mb,Ma)

            if abs(v[a]-v[b]) <= 1e-13
                θab = 0.5*(diag[a]+diag[b])/2.0
            else
                θab = 0.5*(epos[a]-epos[b])/(v[a]-v[b])
            end
            I⁺ += θab*symmetric(Gab+Gba)
        end
    end

    return I⁺
end
=#