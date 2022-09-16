
#Maybe skip AbstractContactSearchAlgorithm and AbstractContactTreatment and just have AbstractSearchAlgorithm
abstract type AbstractContactSearchAlgorithm end
abstract type AbstractContactTreatment end
#abstract type AbstractContactAlgorithm end
#Maybe change name to abstract geometric entity?
abstract type AbstractContactEntity end


struct RigidContactEntity{T} <: AbstractContactEntity
    dofs::Vector{Int}
    segments::Vector{Tuple{Vec{2,T},Vec{2,T}}}
end

#Generic contact entity for all finite elements
struct FaceContactEntity{dim,I,N} <: AbstractContactEntity
    faceinterpolation::I
    dofs::NTuple{N,NTuple{dim,Int}} # N basefunctions
end

function getAABB(a::FaceContactEntity{2,I,N}, x::Vector{T}) where {I,T,N} 
    maxx = maxy = T(-Inf)
    minx = miny = T(Inf)
    for n in 1:N
        k = x[a.dofs[n][1]]
        maxx = k > maxx ? k : maxx
        minx = k < minx ? k : minx

        k = x[a.dofs[n][2]]
        maxy = k > maxy ? k : maxy
        miny = k < miny ? k : miny
    end

    return AABB(Vec{2,T}((minx,miny)), Vec{2,T}((maxx-minx, maxy-miny)))
end

function getAABB(a::FaceContactEntity{3,I,N}, x::Vector{T}) where {T,I,N} 
    maxx = maxy = maxz = T(-Inf)
    minx = miny = minz = T(Inf)
    for n in 1:N
        k = x[a.dofs[n][1]]
        maxx = k > maxx ? k : maxx
        minx = k < minx ? k : minx

        k = x[a.dofs[n][2]]
        maxy = k > maxy ? k : maxy
        miny = k < miny ? k : miny

        k = x[a.dofs[n][3]]
        maxz = k > maxz ? k : maxz
        minz = k < minz ? k : minz
    end

    return AABB(Vec{3,T}((minx,miny,minz)), Vec{3,T}((maxx-minx, maxy-miny, maxz-minz)))
end

function getAABB(a::FaceContactEntity{3,Lagrange{2,RefTetrahedron,1}}, x::Vector{T}) where T
 n1 = x[a.dofs[1]]#assume
 n2 = x[a.dofs[2]]
 n3 = x[a.dofs[3]]
 
 minx, maxx = minimum([n1[1],n2[1],n3[1]]), maximum([n1[1],n2[1],n3[3]])
 miny, maxy = minimum([n1[2],n2[2],n3[2]]), maximum([n1[2],n2[2],n3[2]])
 minz, maxz = minimum([n1[3],n2[3],n3[3]]), maximum([n1[3],n2[3],n3[3]])

 return AABB(Vec{3,T}((minx,miny,minz)), Vec{3,T}((maxx-minx, maxy-miny, maxz-minz)))
end

struct SphericalContactEntity{dim,T} <: AbstractContactEntity
    dofs::Vec{dim,Int}
    radius::T
end

function getAABB(a::SphericalContactEntity{2}, x::Vector{T}) where T 
 x = x[a.dofs]
 r = a.radius
 return AABB(Vec{2,T}((x[1]-r, x[2]-r)), Vec{2,T}((r*2, r*2)))
end

function getAABB(a::SphericalContactEntity{3,T}, x) where T
 x = x[a.dofs]
 r = a.radius
 return AABB(Vec{3,T}((x[1]-r, x[2]-r, x[3]-r)), Vec{3,T}((r*2, r*2, r*2)))
end

struct SegmentContactEntity <: AbstractContactEntity
    dofs1::Vector{Int}
    dofs2::Vector{Int}
end

struct StaticPlaneContactEntity <: AbstractContactEntity
    origin::Vec{2,Float64}
    normal::Vec{2,Float64} #dim,T
end

struct NodeContactEntity{dim} <: AbstractContactEntity
    dofs::NTuple{dim,Int}
    #dofs::Vec{dim,Int}
end

function getAABB(a::NodeContactEntity{dim}, x::AbstractVector{T}) where {dim,T}
 pos = x[a.dofs]#Tuple([x[a.dofs[d]] for d in 1:dim])
 return AABB(pos, zero(Vec{dim,T}))
end

#=
mutable struct Contact_Node2Segment{dim,T} <: AbstractContactSearchAlgorithm
	
	#slaves::Vector{AbstractContactEntity}
	#masters::Vector{AbstractContactEntity}

    #slave_aabb::Vector{AABB{dim,T}}
    #master_aabb::Vector{AABB{dim,T}}

    #ordering::Matrix{Int}  # (antal slavnoder x 5 närmaste mastersegmenten)
    #contact_cache::Vector{Int}

    #slavebuckets::Buckets{dim,T,Int}
    #masterbuckets::Buckets{dim,T,Int}

end
=#
struct ContactOutput{S,M,dim,T}
    slave::S
    master::M
    xi::Vec{2,T} #in 2d, only the first index is used
    penatration::T 
    normal::Vec{dim,T}
    seg_len1::T
    seg_len2::T
    dir1::Vec{dim,T}
    dir2::Vec{dim,T}
end

function search_contact(slave::NodeContactEntity{2}, master::SegmentContactEntity, x)
    dim = 2
    T = Float64

    xs = Vec{dim,T}(Tuple(x[slave.dofs]))

    masterset = zeros(Vec{2, Float64},4)
    x1 = Vec{dim,T}(Tuple(x[master.dofs1]))
    x2 = Vec{dim,T}(Tuple(x[master.dofs2]))

    construct_polygon_from_2d_segment!(masterset, x1, x2)  
    isinside, _segment, _eps, _penatration, _normal, _dir, _seglen = pointInsidePolygon(xs, masterset) 

    if isinside
        return ContactOutput(slave, master, _eps, _penatration, _normal, _seglen, _dir)
    else
        return nothing
    end
end

function search_contact(slave::NodeContactEntity{dim}, master::SphericalContactEntity{dim,T}, x, should_ouput = true) where {dim,T}
    xm = x[master.dofs]
    xs = x[slave.dofs]
    
    dir = (xs-xm)
    dist = norm(dir)
    _normal = dir/dist
    
    _penatration = dist - master.radius 
    
    if _penatration < 0.0
        return ContactOutput(slave, master, Vec{2,T}((T(0.0),T(0.0))), _penatration, _normal, T(0.0), T(0.0), zero(Vec{3,T}), zero(Vec{3,T}))
    else
        return nothing
    end
end

function search_contact(slave::NodeContactEntity{2}, master::FaceContactEntity{2,Lagrange{1,RefCube,1},2}, x::Vector)
    dim = 2
    xs = x[Vec(slave.dofs)]

    x1 = x[Vec(master.dofs[1])]
    x2 = x[Vec(master.dofs[2])]

    _dir = x2 - x1
    _dirlen = norm(_dir)
    _dir /= _dirlen
    
    _normal = [_dir[2], -_dir[1]]
    _normal /= norm(_normal)
    
    d = xs - x1;
    _penatration = dot(d, _normal);
    _eps = dot(d, _dir)/_dirlen;
    
    co = ContactOutput(slave, master, Vec{2}((_eps,0.0)), _penatration, Vec{2}(Tuple(_normal)), _dirlen, 0.0, Vec{2}(Tuple(_dir)), Vec{2}((0.0,0.0)))

    if (_penatration < 0.0 && abs(_penatration) < (_dirlen*0.05)) && (_eps >= 0 && _eps <= 1)
        return true, co
    else
        return false, co
    end
end

function search_contact(slave::NodeContactEntity{3}, master::FaceContactEntity{3,Lagrange{2,RefTetrahedron,1},3}, x::Vector{T}, should_output::Bool = true) where {T}
    dim = 3

    x1 = x[master.dofs[1]]
    
    x2 = x[master.dofs[2]]
    x3 = x[master.dofs[3]]
    
    xs = x[slave.dofs]

    a1 = x2 - x1
    a2 = x3 - x1
    a1_len = norm(a1)
    a2_len = norm(a2)
    a1bar = a1/a1_len
    a2bar = a2/a2_len

    A = [dot(a1,a1) dot(a1,a2); dot(a2,a1) dot(a2,a2)]
    b = [dot(xs-x1, a1), dot(xs-x1, a2)]
    xi = A\b
    xi = Vec{2,T}((xi[1], xi[2]))
    normal = cross(a1bar,a2bar)
    normal /= norm(normal)

    penatration = dot(xs-x1, normal)

    N1 = Ferrite.value(master.faceinterpolation, 1, xi)
    N2 = Ferrite.value(master.faceinterpolation, 2, xi)
    N3 = Ferrite.value(master.faceinterpolation, 3, xi)

    (N1 >=0 && N1 <= 1) || return nothing
    (N2 >=0 && N2 <= 1) || return nothing
    (N3 >=0 && N3 <= 1) || return nothing

    if penatration < 0.0 && abs(penatration) < (a1_len*0.1)
        return ContactOutput(slave, master, xi, penatration, Vec{3,T}(Tuple(normal)), a1_len,a2_len, Vec{3}(Tuple(a1bar)), Vec{3}(Tuple(a2bar)))
    else
        return nothing
    end
end

function search_contact(slave::NodeContactEntity{3}, master::FaceContactEntity{3,Lagrange{2,RefCube,1},4}, x::Vector{T}, ξ::Vec{2,T}=zero(Vec{2,T}), should_outut = true) where {T}
    N = 4
    dim = 3
    @timeit "init search" begin
    nbasefunks = getnbasefunctions(master.faceinterpolation)

    xs = x[slave.dofs]
    xms = [x[master.dofs[i]] for i in 1:nbasefunks]

    #ξ = zero(Vec{2,T})# + Vec{2,T}((0.0,0.0))
    normal = zero(Vec{dim,T})
    a = zeros(Vec{dim,T},2)
    dξ = zeros(T,2)
    xm = zero(Vec{dim,T})
    N = zeros(T, nbasefunks)
    dN = zeros(Tensor{1,2,T}, nbasefunks)
    ddN = zeros(Tensor{2,2,T}, nbasefunks)

    end 

    itr = 0
    @timeit "Find xi" while true
        itr += 1

        @timeit "basefunks" begin
        for i in 1:nbasefunks
            N[i] = Ferrite.value(master.faceinterpolation, i, ξ) 
            #dN[i] = gradient(ξ -> Ferrite.value(master.faceinterpolation, i, ξ), ξ)
            #ddN[i] = hessian(ξ -> Ferrite.value(master.faceinterpolation, i, ξ), ξ)
        end

        dN[1] = 0.25*Vec(-(1-ξ[2]), -(1-ξ[1]))
        dN[2] = 0.25*Vec((1-ξ[2]), -(1+ξ[1]))
        dN[3] = 0.25*Vec((1+ξ[2]), (1+ξ[1]))
        dN[4] = 0.25*Vec(-(1+ξ[2]), (1-ξ[1]))

        ddN[1] = SymmetricTensor{2,2,T}((0.0,0.25,0.0))
        ddN[2] = SymmetricTensor{2,2,T}((0.0,-0.25,0.0))
        ddN[3] = SymmetricTensor{2,2,T}((0.0,0.25,0.0))
        ddN[4] = SymmetricTensor{2,2,T}((0.0,-0.25,0.0)) 

        end

        @timeit "solve" begin
        
        fill!(dξ,0.0)
        rhs = 0.0
        lhs = 0.0
        for α in 1:2
            lhs = 0.0
            rhs = 0.0
            dxm = zero(Vec{dim,T})
            xm = zero(Vec{dim,T})

            for l in 1:nbasefunks
                xm += N[l]*xms[l]
            end

            
            normal = xs - xm

            
            for n in 1:nbasefunks
                dxm += dN[n][α]*xms[n]
            end
            a[α] = dxm
            for β in 1:2
                for i in 1:nbasefunks
                    for j in 1:nbasefunks
                        lhs += dN[i][α]*dN[j][β]*(xms[i]⋅xms[j])
                        lhs -= ddN[j][α,β]*(xms[j]⋅normal)
                    end
                end
            end
            rhs = (normal)⋅dxm
            dξ[α] = rhs/lhs

        end
        end

        if itr > 10
            @timeit "return" co = ContactOutput(slave, master, ξ, -1.0, zero(Vec{3,T}), -1.0, -1.0, zero(Vec{3,T}), zero(Vec{3,T}))
            return false, co
        end
        if norm(dξ) < 0.0001
            break
        end

        dξ2 = Vec{2,T}((dξ[1],dξ[2]))
        ξ += dξ2
    end
    
    normal = cross(a[1],a[2])
    normal /= norm(normal)
    
    penatration = dot(xs-xm, normal)

    normal *= sign(penatration)
    
###
    
    if abs(penatration) > 2.0
        return false, ContactOutput(slave, master, ξ, -1.0, zero(Vec{3,T}), -1.0, -1.0, zero(Vec{3,T}), zero(Vec{3,T}))
    end
    penatration = -(2 - abs(penatration))
    #Check if inside bounds
    for i in 1:nbasefunks
        if N[i] < 0.0 || N[i] > 1.0
            return false, ContactOutput(slave, master, ξ, -1.0, zero(Vec{3,T}), -1.0, -1.0, zero(Vec{3,T}), zero(Vec{3,T})) 
        end
    end

    return true, ContactOutput(slave, master, ξ, penatration, normal, -1.0, -1.0, zero(Vec{3,T}), zero(Vec{3,T}))

###

    #=if abs(penatration) > 2.0
        return nothing
    end
    penatration = -(2 - abs(penatration))
    #Check if inside bounds
    for i in 1:nbasefunks
        if N[i] < 0.0 || N[i] > 1.0
            return nothing
        end
    end

    return ContactOutput(slave, master, ξ, penatration, normal, -1.0, -1.0, zero(Vec{3,T}), zero(Vec{3,T}))
    =#
end

#
# Contact treatment, probable move to other file later
#

struct PenaltyBasedContactWithoutFriction{T} <: AbstractContactTreatment
    penalty::T
    #friction::T
end

function handle_contact!(co::ContactOutput{S,M,dim,T}, ct::PenaltyBasedContactWithoutFriction, f, K, x) where {S<:NodeContactEntity, M<:SphericalContactEntity, dim, T}
    #f[co.slave.dofs]  = co.normal*co.penatration*ct.penalty

    xm = x[co.master.dofs]
    xs = x[co.slave.dofs]
    

    #dofs = vcat(co.master.dofs[1:dim], co.slave.dofs[1:dim])
    F = co.penatration
    G = co.normal
    #H = ForwardDiff.hessian((x) -> force_on_rigidbody(x), x[dofs])
    f[co.master.dofs] += ct.penalty*F*G
    f[co.slave.dofs] -= ct.penalty*F*G
    #K[dofs,dofs] += (ct.penalty*G*G' + ct.penalty*H*F)
end

function handle_contact!(co::ContactOutput{S,M,dim,T}, ct::PenaltyBasedContactWithoutFriction, f, K, x) where {S<:NodeContactEntity, M<:StaticPlaneContactEntity, dim, T}
    
    f[co.slave.dofs] -= ct.penalty*co.penatration*co.master.normal
end

function handle_contact!(co::ContactOutput{S,FaceContactEntity{2,Lagrange{1,RefCube,1},2},dim,T}, ct::PenaltyBasedContactWithoutFriction, f::Vector{T}) where {S<:NodeContactEntity, dim, T}
    
    N = co.normal
    l = co.seg_len1
    xi = co.xi[1]
    gns = co.penatration

    #Ns = vcat(N, -(1-xi)*N, -xi*N)
    #N0s = vcat([0.0,0.0], -N, N)

    f[Vec(co.slave.dofs)] .-= ct.penalty*gns*N
    f[Vec(co.master.dofs[1])] .+= ct.penalty*gns*(1-xi)*N
    f[Vec(co.master.dofs[2])] .+= ct.penalty*gns*xi*N

    #co_dofs = vcat(co.slave.dofs, co.master.dofs)

    #f[co_dofs] += -ct.penalty*gns*Ns
    
end

function handle_contact!(co::ContactOutput{S,FaceContactEntity{3,Lagrange{2,RefTetrahedron,1},3},dim,T}, ct::PenaltyBasedContactWithoutFriction, f::Vector{T}, K, x) where {S<:NodeContactEntity, dim, T}
    xi = co.xi
    n = co.normal
    penatration = co.penatration

    N1 = Ferrite.value(co.master.faceinterpolation, 1, xi)
    N2 = Ferrite.value(co.master.faceinterpolation, 2, xi)
    N3 = Ferrite.value(co.master.faceinterpolation, 3, xi)

    #println("Slave contact master, pen: $penatration, slavedof: $(co.slave.dofs)")

    f[co.slave.dofs]       -= ct.penalty*penatration*n
    f[co.master.dofs[1]] += ct.penalty*penatration*N1*n
    f[co.master.dofs[2]] += ct.penalty*penatration*N2*n
    f[co.master.dofs[3]] += ct.penalty*penatration*N3*n
    
end

function handle_contact!(co::ContactOutput{S,FaceContactEntity{3,Lagrange{2,RefCube,1},4},dim,T}, ct::PenaltyBasedContactWithoutFriction, f::Vector{T}, K, x) where {S<:NodeContactEntity, dim, T}
    n = co.normal
    penatration = co.penatration

    #The basefunction values should be stored in ContactOutput
    N = Ferrite.value(co.master.faceinterpolation, co.xi)

    f[co.slave.dofs]       -= ct.penalty*penatration*n
    for i in 1:getnbasefunctions(co.master.faceinterpolation)
        f[co.master.dofs[i]] += ct.penalty*penatration*N[i]*n
    end
    
end

struct ContactHandler
    search_algo::AbstractContactSearchAlgorithm #AbstractSearchAlgo
    contact_treatment::AbstractContactTreatment
end

function contact!(search_algo, state, globaldata)
    x = state.x0 .+ state.d
    
    @timeit "Updating contact" update_contact!(search_algo, x)
    @timeit "Searching contact" contactoutputs, ncontacts = search1!(search_algo, x)

    contact_treatment = PenaltyBasedContactWithoutFriction(2e5)
    for i in 1:ncontacts#contactoutputs
        co = contactoutputs[i]
        
        @timeit "Handling contact" handle_contact!(co, contact_treatment, state.system_arrays.fᵉ)
    end
end

#
# Utilitiy
# 
#=
function pointInsidePolygon(point, vertices::Vector{Vec{2, Float64}})# where dim

	dim = 2
	nverts = length(vertices)

	inside = false
	min_distance = Inf;
	
	_penatration = Inf; 
	_eps = 0;
	_seglen = 0;
	_segment = 0;
	_normal = zeros(Float64,dim)
	_dir = zeros(Float64,dim)
	
	for iv in 1:nverts
		ivert = vertices[iv]
		jvert = vertices[iv+1 > nverts ? 1 : iv+1]
		
		dir = ivert - jvert;
		dirlen = norm(dir);
		dir = dir/dirlen;
		
		normal = Vec{dim, Float64}((-dir[2], dir[1]));
		normal /= norm(normal);
		
		d = point - jvert;
		
		ipenatration = dot(d, normal);
		ieps = dot(d, dir);
		
		if(ipenatration[1] > 0.0)
			inside = false;
			return inside, nothing, nothing, nothing, nothing, nothing, nothing
		end
		if abs(ipenatration[1]) < abs(_penatration)
			if norm(d) < min_distance
				min_distance = norm(d);
				_penatration = ipenatration[1];
				_eps = ieps;
				_normal = normal;
				_segment = iv
				_seglen = dirlen;
				inside = true;
				_dir = dir;
			end
		end
	end

	return inside, _segment, _eps, _penatration, _normal, _dir, _seglen
end

function construct_polygon_from_2d_segment!(polygon, x1, x2, tt = 1)

    #The mastersegment will be a line with two nodes 
    #tt = 0.1; #Some sort of thickness
    #x1 = x[mastersegment[1]].x
    #x2 = x[mastersegment[2]].x
    
    dir = x1 - x2;
    dirlen = norm(dir);
    dir = dir/dirlen;
    normal = Vec{2,Float64}((-dir[2], dir[1]));
    normal /= norm(normal);

    polygon[1] = (x1 + normal*tt)
    polygon[2] = (x2 + normal*tt)
    polygon[3] = (x2 - normal*tt)
    polygon[4] = (x1 - normal*tt)

end
=#