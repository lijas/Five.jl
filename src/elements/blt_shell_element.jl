
"""
BelytschkoLinTsayShellElement
http://roadsafellc.com/NCHRP22-24/Literature/Papers/Explicit%20Algorithms%20for%20the%20nonlinear%20dynamics%20of%20shells.pdf
Explicit Algorithms for the nonlinear dynamics of shells, Ted belytschko and jerr I lin
"""

struct BelytschkoLinTsayShellElementInfo{T}
    thickness::T
    shear_factor::T
    nip::Int
end

struct BelytschkoLinTsayShellElement{T} <: AbstractElement
    
    info::BelytschkoLinTsayShellElementInfo{T}

    fields::Vector{Field}
    celltype::Type{<:Cell}

    cv::CellScalarValues{2,T,RefCube}
    ip::QuadratureRule{1,RefCube,T} #Quadpoints through thickness

    v̂ⁿ::Vector{Vec{3,T}}
    θ̂ⁿ::Vector{Vec{3,T}}
    B::Matrix{T}
    N::Vector{T}

    fᴿ::Matrix{T}
    mᴿ::Matrix{T}

    #R::Tensor{2,3,T,9} #rotation matrix
    A::Ref{T} #Area
    x̂_2d::Vector{Vec{2,T}}
end

mutable struct BelytschkoLinTsayShellElementState{T} <: AbstractElementState
    ue::Vector{T}
    Qᴹ::Vector{T}
    Qᵇ::Vector{T}
end

get_elementstate_type(e::BelytschkoLinTsayShellElement) = BelytschkoLinTsayShellElementState
function BelytschkoLinTsayShellElementState(e::BelytschkoLinTsayShellElement{T}, coords) where T
    return BelytschkoLinTsayShellElementState{T}(zeros(T,24), zeros(T,2), zeros(T,3))
end

JuAFEM.getnquadpoints(e::BelytschkoLinTsayShellElement) = 5
JuAFEM.ndofs(e::BelytschkoLinTsayShellElement) = 24
has_constant_massmatrix(::BelytschkoLinTsayShellElement) = true

function BelytschkoLinTsayShellElement{T}(;thickness, shear_factor=5/6, nip=5) where {T}
    
    ip = Lagrange{2,RefCube,1}()
    qr = QuadratureRule{2,RefCube}(1)
    

    fields = [Field(:u, Lagrange{2,RefCube,1}(), 3), 
              Field(:θ, Lagrange{2,RefCube,1}(), 3)]
    cv = CellScalarValues(qr,ip) 
    
    info = BelytschkoLinTsayShellElementInfo{T}(thickness, shear_factor, nip)

    return BelytschkoLinTsayShellElement{T}(info, fields, Cell{3,4,2}, cv, QuadratureRule{1,RefCube}(5), zeros(Vec{3,T},4), zeros(Vec{3,T},4), zeros(T,4,2), zeros(T,4), zeros(T,3,3),zeros(T,3,3), Ref{T}(0.0), zeros(Vec{2,T},4))
end

function calculate_minimum_timestep(element::BelytschkoLinTsayShellElement{T}, material::AbstractMaterial, cell::CellIterator, ue::Vector, due::Vector) where {T}
    error("Not implemented")
end

function integrate_forcevector!(element::BelytschkoLinTsayShellElement{T}, elementstate::BelytschkoLinTsayShellElementState, material::AbstractMaterial, materialstate::Vector{<:AbstractMaterialState}, fe::Vector, cell, Δue::Vector, ue::Vector, due::Vector) where {T}
    
    nnodes = 4
    ndofs = 24

    κ = element.info.shear_factor
    θ̂ⁿ = element.θ̂ⁿ
    v̂ⁿ = element.v̂ⁿ
    B  = element.B
    N = element.N
    fᴿ = element.fᴿ
    mᴿ = element.mᴿ
    ip = element.ip
    #deformation since last timestep
    #Δue = ue - elementstate.ue
    #elementstate.ue = copy(ue)

    u_range = 1:(nnodes*3)
    θ_range = u_range .+ nnodes*3

    u = reinterpret(Vec{3,T}, ue[u_range])
    #v = reinterpret(Vec{3,T}, due[u_range])
    #θ = reinterpret(Vec{3,T}, due[θ_range])

    x = cell .+ u
    Δx = reinterpret(Vec{3,T}, Δue[u_range])
    Δθ = reinterpret(Vec{3,T}, Δue[θ_range])
    
    
    R = _reinit_blt!(element, x, Δx, Δθ)
    #R = element.R
    A = element.A[]
    #ip = QuadratureRule{1,RefCube}(5)

    fill!(fᴿ, 0.0)# = zeros(T, 3, 3)
    fill!(mᴿ, 0.0)# = zeros(T, 3, 3)

    # HOURGLASS
    #=
    γ= zeros(T,4)
    h = [1,-1,1,-1]
    for i in 1:4

        #Hourglass vector γ
        _hxb = 0
        for j in 1:4
            for d in 1:2
                _hxb += h[j].*element.x̂_2d[j][d]*element.B[i,d]
            end
        end
        γ[i] = h[i] - _hxb
    end
    
    
    qᵇ = [sum([γ[i]*θ̂ⁿ[i][1] for i in 1:4]),
          sum([γ[i]*θ̂ⁿ[i][2] for i in 1:4]),
          sum([γ[i]*v̂ⁿ[i][3] for i in 1:4])]

    qᴹ = [sum([γ[i]*v̂ⁿ[i][1] for i in 1:4]),
          sum([γ[i]*v̂ⁿ[i][2] for i in 1:4])]
    =#
    vxx = vyy = vxy = vxz = vyz = vzy = vyx = vzx = 0.0
    for iqp in 1:length(ip.weights)

        zeta = ip.points[iqp][1]*element.info.thickness
        dz   = ip.weights[iqp]*element.info.thickness*0.5
        
        #=vxx = dv̂ᵐdx̂[1,1] + zeta*dθ̂dx̂[2,1]
        vyy = dv̂ᵐdx̂[2,2] - zeta*dθ̂dx̂[1,2]
        vxy = dv̂ᵐdx̂[1,2] + zeta*dθ̂dx̂[2,2]
        vyx = dv̂ᵐdx̂[2,1] - zeta*dθ̂dx̂[1,1]

        vyz = -θ̂[1]
        vzy = dv̂ᵐdx̂[3,2]

        vxz = +θ̂[2]
        vzx = dv̂ᵐdx̂[3,1]=#
        vxx = vyy = vxy = vxz = vyz = vzy = vyx = vzx = 0.0
        @inbounds for i in 1:4
            vxx += B[i,1]*v̂ⁿ[i][1] + zeta*B[i,1]*θ̂ⁿ[i][2]
            vyy += B[i,2]*v̂ⁿ[i][2] - zeta*B[i,2]*θ̂ⁿ[i][1]
            vxy += B[i,2]*v̂ⁿ[i][1] + zeta*B[i,2]*θ̂ⁿ[i][2]
            vyx += B[i,1]*v̂ⁿ[i][2] - zeta*B[i,1]*θ̂ⁿ[i][1]
            
            vyz += -N[i]*θ̂ⁿ[i][1]
            vzy += B[i,2]*v̂ⁿ[i][3]

            vxz += N[i]*θ̂ⁿ[i][2]
            vzx += B[i,1]*v̂ⁿ[i][3]
        end
        velocity_gradient = Tensor{2,3,T}((vxx,vxy,vxz,vyx,vyy,vyz,vzx,vzy,0.0))
        
        #iterate for plane stress
        TOL = 0.001
        normr = 1
        local sigma, new_matstate
        for i in 1:3
            sigma, C, new_matstate = constitutive_driver(material, velocity_gradient, materialstate[iqp])

            dS3333 = C[3,3,3,3]
            res = -sigma[3,3]
            normr = norm(res)

            dE = inv(dS3333)*res
            velocity_gradient += Tensor{2,3,T}((0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,dE))
        end
        materialstate[iqp] = new_matstate

        #@show dz, zeta
        @inbounds for i in 1:3
             for j in 1:3
                fᴿ[i,j] += sigma[i,j] * dz
                mᴿ[i,j] -= zeta*sigma[i,j] * dz
            end
        end
    end

    
    #Hourglass stress
    #=
    C1,C2,C3 = (0.0,0.0,0.0)
    _BB = 0
    rθ = 0.01
    rM = 0.01
    rW = 0.01
    for i in 1:4
        for d in 1:2
            _BB += B[i,d]
        end
    end

    E = material.E
    G = material.E/2/(material.nu+1)

    C1 = 1/192*rθ*(A*E*elementinfo.thickness^3)*_BB
    C2 = 1/12 *rW*(κ*G*elementinfo.thickness^3)*_BB
    C3 = 1/8  *rM*(A*G*elementinfo.thickness  )*_BB

    Qᵇ = vcat(C1*qᵇ[1:2], 
              C2*qᵇ[3]) + elementstate.Qᵇ
    Qᴹ = C3*qᴹ + elementstate.Qᴹ 

    elementstate.Qᴹ = Qᴹ
    elementstate.Qᵇ = Qᵇ
    =#

    counter = 1
   @inbounds for i in 1:nnodes
        f̂ = Vec( A*(B[i,1]*fᴿ[1,1] + B[i,2]*fᴿ[1,2])  ,#   +   γ[i]*Qᴹ[1],
                 A*(B[i,2]*fᴿ[2,2] + B[i,1]*fᴿ[1,2])      ,#  +   γ[i]*Qᴹ[2],
                 A*κ*(B[i,1]*fᴿ[1,3] + B[i,2]*fᴿ[2,3])    )#  +   γ[i]*Qᵇ[3])

        m̂ = Vec((
                A*(B[i,2]*mᴿ[2,2] + B[i,1]*mᴿ[1,2] - κ*0.25*fᴿ[2,3])   ,#  +   γ[i]*Qᵇ[1],
                A*(-B[i,1]*mᴿ[1,1] - B[i,2]*mᴿ[1,2] + κ*0.25*fᴿ[1,3])   ,# +   γ[i]*Qᵇ[2],
                0.0))

        fi::Vec{3,T} = R'⋅f̂ 
        mi::Vec{3,T} = R'⋅m̂


        @inbounds for d in 1:3 #dim
            fe[counter] = fi[d]
            fe[counter+12] = mi[d]
            counter += 1
        end
    end

end

function integrate_massmatrix!(element::BelytschkoLinTsayShellElement{T}, elstate::BelytschkoLinTsayShellElementState, material::AbstractMaterial, cell::CellIterator, me::Matrix, ue::AbstractVector, due::AbstractVector) where {T}

    nnodes = 4
    cv = element.cv

    x = cell.coords

    x̂_2d = element.x̂_2d

    R = _reinit_blt!(element, x, zeros(Vec{3,T}), zeros(Vec{3,T})) 
    A = element.A
    reinit!(cv, x̂_2d)

    area = 0

    _me = zeros(T, nnodes, nnodes)
    for qp in 1:getnquadpoints(cv)
        dA = getdetJdV(cv, qp)
        
        for i in 1:nnodes
            Ni = shape_value(cv, qp, i)
            for j in 1:nnodes
                Nj = shape_value(cv, qp, j)
                _me[i,j] += density(material)*Ni⋅Nj*dA
                area += dA
            end
        end
    end

    _me .*= elinfo.thickness
    volume = area*elinfo.thickness
    
    _me_lumped = zeros(T, nnodes, nnodes)
    for i in 1:nnodes
        for j in 1:nnodes
            _me_lumped[i,i] += _me[i,j]
        end
    end

    a1 = (1/12)*elinfo.thickness^2
    a2 = volume/8/elinfo.thickness
    alpha = max(a1,a2)
    #alpha = 100
    c=1
    for i in 1:nnodes
        for d in 1:3 #dim
            me[c,c] = _me_lumped[i,i]
            me[c+nnodes*3,c+nnodes*3] = alpha*_me_lumped[i,i]
            c+=1
        end
    end

end

function integrate_forcevector_and_stiffnessmatrix!(element::BelytschkoLinTsayShellElement{T}, elementstate::BelytschkoLinTsayShellElementState{T}, material::AbstractMaterial, materialstate::Vector{AbstractMaterialState}, fe::Vector, ke::AbstractMatrix, cell, ue::Vector, due::Vector) where {T}
   

end

function bodyforce!(element::BelytschkoLinTsayShellElement{T}, elstate::BelytschkoLinTsayShellElementState, material::AbstractMaterial, cell::CellIterator, fe::Vector, forcevec::Vector) where T
    #@assert length(forcevec) == dim

    x̂_2d = element.x̂_2d
     
    _reinit_blt!(element, cell.coords, zeros(Vec{3,T}), zeros(Vec{3,T}))

    cv = element.cv
    reinit!(cv, x̂_2d)
    
    area = 0
    for qp in 1:getnquadpoints(cv)
        area += getdetJdV(cv, qp)
    end
    
    volume = area*element.info.thickness

    c = 0
    for i in 1:4
        force = forcevec*density(material)*volume/4
        for j in 1:3
            c += 1
            fe[c] = force[j]
        end 
    end

end

@inline function _reinit_blt!(element::BelytschkoLinTsayShellElement, x::Vector{Vec{3,T}}, Δx, Δθ) where T

    nnodes = 4

    r21 = x[2] - x[1]
    r31 = x[3] - x[1]
    r42 = x[4] - x[2]
    s3 = cross(r31, r42)
    e3 = s3/norm(s3)

    s1 = r21 - dot(r21, e3)*e3
    e1 = s1/norm(s1)

    e2 = cross(e3,e1)

    R = Tensor{2,3,T}((e1[1],e2[1],e3[1], 
                       e1[2],e2[2],e3[2],
                       e1[3],e2[3],e3[3]))

    #x̂_3d = zeros(Vec{3,T}, nnodes)
    #@inbounds for i in 1:nnodes
    #    x̂_3d[i] = R⋅x[i]
    #end
    
    #Ignore z-component
    x̂_2d = element.x̂_2d
    @inbounds for i in 1:nnodes
        x̂_3d = R⋅x[i]
        x̂_2d[i] = Vec{2,T}( (x̂_3d[1], x̂_3d[2]) )
    end 

    
    A = 0.5*((x̂_2d[3]-x̂_2d[1])[1]*(x̂_2d[4]-x̂_2d[2])[2] - (x̂_2d[2]-x̂_2d[4])[1]*(x̂_2d[1]-x̂_2d[3])[2])
    element.A[] = A
    
    #B2 = 1/2/A[]*[x̂_2d[2][2]-x̂_2d[4][2] x̂_2d[3][2]-x̂_2d[1][2] x̂_2d[4][2]-x̂_2d[2][2] x̂_2d[1][2]-x̂_2d[3][2]; 
    #                    x̂_2d[4][1]-x̂_2d[2][1] x̂_2d[1][1]-x̂_2d[3][1] x̂_2d[2][1]-x̂_2d[4][1] x̂_2d[3][1]-x̂_2d[1][1]]' 
    
    
    element.B[1,1] = x̂_2d[2][2]-x̂_2d[4][2] 
    element.B[2,1] = x̂_2d[3][2]-x̂_2d[1][2] 
    element.B[3,1] = x̂_2d[4][2]-x̂_2d[2][2] 
    element.B[4,1] = x̂_2d[1][2]-x̂_2d[3][2] 
    element.B[1,2] = x̂_2d[4][1]-x̂_2d[2][1] 
    element.B[2,2] = x̂_2d[1][1]-x̂_2d[3][1] 
    element.B[3,2] = x̂_2d[2][1]-x̂_2d[4][1] 
    element.B[4,2] = x̂_2d[3][1]-x̂_2d[1][1] 

    element.B .*= 1/2/A
    
    element.N .= (0.25,0.25,0.25,0.25)

    @inbounds for i in 1:nnodes
        element.v̂ⁿ[i] = R⋅Δx[i]
        element.θ̂ⁿ[i] = R⋅Δθ[i]
    end

    return R

end

#In 2d, the face interpolation type is a "1d line"
#getfaceinterpolation(element::BelytschkoLinTsayShellElement{T}) where {T} = Lagrange{1, RefCube, 1}()
#getfaceinterpolation(element::BelytschkoLinTsayShellElement{T}) where {T} = Lagrange{2, shape, 1}()


