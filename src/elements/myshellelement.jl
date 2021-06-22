
"""
ShellElement


"""

struct ShellElement{T} <: AbstractElement
    ndofs::Int
    fields::Vector{Field}
    celltype::Type{<:Cell}

    #Probably only need one of theese. Fix later...
    cv::CellScalarValues{3,T,RefCube}
    cv2::CellVectorValues{3,T,RefCube,9}
end

struct ShellElementInfo{T}
    thickness::T
    shear_factor::T
    nip::Int
end

mutable struct ShellElementState{T} <: AbstractElementState
    p::Vector{Vec{3,T}}
    θ::Vector{Vector{T}}
    X::Vector{Vec{3,T}} #original coords of solid element
end

get_elementinfo_type(e::ShellElement) = ShellElementInfo
get_elementstate_type(e::ShellElement) = ShellElementState
function ShellElementState(e::ShellElement{T}, info::ShellElementInfo, coords) where T
    a1 = coords[2] - coords[1]
    a2 = coords[4] - coords[1]
    a1 /= norm(a1)
    a2 /= norm(a2)
    n = cross(a1,a2)
    n /= norm(n)
    p = [deepcopy(n) for i in 1:4]
    θ = [zeros(T,3) for i in 1:4]#zeros(Vec{3,T},4)
    
    nnodes_shell = 4
    X = zeros(Vec{3,T},nnodes_shell*2)
    for i in 1:nnodes_shell
        X[i] = coords[i] - info.thickness/2*p[i]
        X[i+nnodes_shell] = coords[i] + info.thickness/2*p[i]    
    end

    return ShellElementState(p,θ,X)

end
Ferrite.getnquadpoints(e::ShellElement) = Ferrite.getnquadpoints(e.cv)
Ferrite.ndofs(e::ShellElement) = e.ndofs

function ShellElement{T}() where {T}
    
    ip = Lagrange{3, RefCube, 1}()
    qr = QuadratureRule{3, RefCube}(2)
    cv_solid = CellScalarValues(qr,ip)

    fields = [Field(:u, Lagrange{2, RefCube, 1}(), 3), 
              Field(:θ, Lagrange{2, RefCube, 1}(), 3)]
    cv2 = CellVectorValues(qr,ip) 

    return ShellElement{T}(6*4, fields, Cell{3,4,2}, cv_solid, cv2)
end

function calculate_minimum_timestep(element::ShellElement{T}, info::ShellElementInfo, material::AbstractMaterial, cell::CellIterator, ue::Vector, due::Vector) where {T}
    error("Not implemented")
end

function integrate_forcevector!(element::ShellElement{T}, elementstate::ShellElementState, elementinfo::ShellElementInfo, material::AbstractMaterial, materialstate::Vector{<:AbstractMaterialState}, fe::Vector, cell, ue::Vector, due::Vector) where {T}

    cv = element.cv2

    nnodes_shell = 4
    h = elementinfo.thickness

    old_dirctors = elementstate.p
    
    directors = zeros(Vec{3,T}, nnodes_shell)
    solid_coords = zeros(Vec{3,T}, nnodes_shell*2)

    un = reinterpret(Vec{3,T}, ue[1:(nnodes_shell*3)])
    θn = reinterpret(Vec{3,T}, ue[(nnodes_shell*3+1):end])
    #un = reinterpret(NTuple{6,T}, ue)
    xᵐ = zeros(Vec{3,T},nnodes_shell)

    scew(v) = Tensor{2,3,T}((0.0, v[3],-v[2], -v[3], 0.0, v[1], v[2], -v[1], 0.0))

    e = (Vec{3,T}((1.0,0.0,0.0)),Vec{3,T}((0.0,1.0,0.0)),Vec{3,T}((0.0,0.0,1.0)))

    for i in 1:nnodes_shell
        #xᵐ[i] = cell[i] + Vec{3,T}(un[i][1:3])
        #dθ = collect(un[i][4:6]) - elementstate.θ[i]
        xᵐ[i] = cell[i] + un[i]
        dθ = θn[i] - elementstate.θ[i]
        R = inv(one(Tensor{2,3,T}) - 0.5*scew(dθ))⋅(one(Tensor{2,3,T}) + 0.5*scew(dθ))
        #R = one(Tensor{2,3,T}) + 1/(1 + 0.25*tt^2)*(scew(dθ) + 0.5*scew(dθ)⋅scew(dθ))
        #θ = norm(dθ)
        #v = dθ/θ
        #v = reinterpret(Vec{3,T}, v)[1]

        #q0, q = cos(θ/2), v*sin(θ/2)
        #R = 2*(q0^2 - 0.5)*one(Tensor{2,3,T}) + 2*q0*scew(q) + 2*q⊗q

        #@show R
        directors[i] = R⋅old_dirctors[i]
        elementstate.p[i] = directors[i]
        elementstate.θ[i] += dθ
    end
    
    #@show directors[2]
    for i in 1:nnodes_shell
        solid_coords[i+nnodes_shell] = xᵐ[i] + h/2*directors[i]
        solid_coords[i] = xᵐ[i] - h/2*directors[i]
    end

    #δE = zeros(SymmetricTensor{2, 3, T, 6}, (nnodes_shell*2)*3)
 
    fs = zeros(Vec{3,T}, nnodes_shell*2)
    fs2 = zeros(T, nnodes_shell*2*3)
    reinit!(cv, elementstate.X)

    ndofs = 8*3
    δE = zeros(SymmetricTensor{2, 3, T, 6}, ndofs)

    #@show solid_coords
    for qp in 1:getnquadpoints(cv)

        #=
        g1 = zero(Vec{3,T})
        g2 = zero(Vec{3,T})
        for i in 1:(nnodes_shell*2)
            G = shape_gradient(cv, qp, solid_coords)
            g1 += cv.dNdξ[i,qp][1]*solid_coords[i]
            g2 += cv.dNdξ[i,qp][2]*solid_coords[i]
        end

        ez = cross(g1,g2)/norm(cross(g1,g2))
        a = g1/norm(g1) + g2/norm(g2)
        b = cross(ez,a)
        ex = (a-b)/norm(a-b)
        ey = (a+b)/norm(a+b)
    

        R = Tensor{2,3,T}( (e[1]⋅ex, e[1]⋅ey, e[1]⋅ez, 
                            e[2]⋅ex, e[2]⋅ey, e[2]⋅ez,
                            e[3]⋅ex, e[3]⋅ey, e[3]⋅ez) )
        =#

        ∇u = function_gradient(cv, qp, solid_coords .- elementstate.X)
        dΩ = getdetJdV(cv, qp)

        F = one(∇u) + ∇u
        E = symmetric(1/2 * (F' ⋅ F - one(F)))

        TOL = 0.001
        normr = 1
        local S, new_matstate
        while normr > TOL
            S, ∂S∂E, new_matstate = constitutive_driver(material, E, materialstate[qp])

            dS3333 = ∂S∂E[3,3,3,3]
            res = -S[3,3]
            normr = norm(res)

            dE = inv(dS3333)*res
            E += Tensor{2,3,T}((0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,dE))
        end
        materialstate[qp] = new_matstate

        #P = S⋅F'
        #σ = inv(J)*F⋅S⋅F'
        #σ = R'*σ*R

        #=for inode in 1:(nnodes_shell*2)
            δFi = shape_gradient(cv, qp, inode)
            fs[inode] += (P ⋅ δFi) * dΩ
        end=#
        for i in 1:ndofs
            δu = shape_value(cv, qp, i)
            δFi = shape_gradient(cv, qp, i)
            δE[i] = symmetric(1/2*(δFi'⋅F + F'⋅δFi))
            fs2[i] += (δE[i] ⊡ S) * dΩ
        end

    end
    fs = reinterpret(Vec{3,T},fs2)
    
    for i in 1:nnodes_shell
        pi = 0.5*h*directors[i]
        #T_ = Tensor{3,2,T}((0,pi[1],pi[3],pi[1]))
        c = (i-1)*3
        f_trans = +fs[i]+fs[i+4]
        f_bend = -scew(pi)⋅fs[i] + scew(pi)⋅fs[i+4]
        #@show f_bend, f_trans
        for j in 1:3
            fe[c+j] = f_trans[j]
            fe[c+j+nnodes_shell*3] = f_bend[j]

        end
    end


end

function integrate_massmatrix!(element::ShellElement{T}, elstate::ShellElementState, elinfo::ShellElementInfo, material::AbstractMaterial, cell::CellIterator, me::Matrix, ue::AbstractVector, due::AbstractVector) where {T}

    cv = element.cv

    reinit!(cv, elstate.X)
    nnodes_shell = 4
    ndofs = nnodes_shell*2*3
    
    volume = 0
    node_mass_solid = zeros(T,nnodes_shell*2)
    for qp in 1:getnquadpoints(cv)
        dV = getdetJdV(cv, qp)
        
        c = 0
        for i in 1:nnodes_shell*2
            c += 1
            Ni = shape_value(cv, qp, i)
            for j in 1:nnodes_shell*2
                Nj = shape_value(cv, qp, j)
                node_mass_solid[c] += density(material)*Ni⋅Nj*dV
                volume +=dV
            end
        end
    end
    
    node_mass_shell = zeros(T,nnodes_shell)
    for i in 1:nnodes_shell
        node_mass_shell[i] = node_mass_solid[i] + node_mass_solid[i+nnodes_shell]
    end


    a1 = (1/12)*elinfo.thickness^2
    a2 = volume/8/elinfo.thickness
    alpha = max(a1,a2)
    #alpha = 100
    c=1
    for i in 1:nnodes_shell
        for d in 1:3
            me[c,c] = node_mass_shell[i]
            me[c+nnodes_shell*3,c+nnodes_shell*3] = alpha*node_mass_shell[i]
            c+=1
        end
    end

    @show alpha 
    error("hej")


end

function integrate_forcevector_and_stiffnessmatrix!(element::ShellElement{T}, elementstate::ShellElementState{T}, elementinfo::ShellElementInfo, material::AbstractMaterial, materialstate::Vector{AbstractMaterialState}, fe::Vector, ke::AbstractMatrix, cell, ue::Vector, due::Vector) where {T}
   

end

function bodyforce!(element::ShellElement{T}, elstate::ShellElementState, elementinfo::ShellElementInfo, material::AbstractMaterial, cell::CellIterator, fe::Vector, forcevec::Vector) where T
    #@assert length(forcevec) == dim
    cv = element.cv
    reinit!(cv, elstate.X)
    
    volume = 0
    for qp in 1:getnquadpoints(cv)
        volume += getdetJdV(cv, qp)
    end
    c = 0
    for i in 1:4
        force = forcevec*density(material)*volume/4
        for j in 1:3
            c += 1
            fe[c] = force[j]
        end 
    end

end

#In 2d, the face interpolation type is a "1d line"
#getfaceinterpolation(element::ShellElement{T}) where {T} = Lagrange{1, RefCube, 1}()
#getfaceinterpolation(element::ShellElement{T}) where {T} = Lagrange{2, shape, 1}()


