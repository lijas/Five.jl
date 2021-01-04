
struct AABB{dim,T}
    cornerpos::Vec{dim,T} #change to Vec{dim,T}
    sidelength::Vec{dim,T}
end

#=function AABB(vars...)
    n = length(vars)
    @assert n%2 == 0
    return AABB(Tuple(vars[1:n/2]), Tuple(vars[(n/2):n]))
end=#

struct Bucket{dim,T<:AbstractFloat,obj}
	bounds::AABB{dim,T}
	objects::Vector{obj}
end

struct Buckets{dim,T,obj}
    objects::Vector{Vector{obj}}
    #buckets::Vector{Bucket{dim,T,obj}}
    worldbounds::AABB{dim,T}
    nbuckets::NTuple{dim,Int}
end

function Buckets(::Type{obj}, worldbounds::AABB{dim,T}, nbuckets::NTuple{dim,Int}) where {dim,T,obj}
    ntotalbuckets = prod(nbuckets)
    buckets = Vector{Vector{obj}}()
    for i in 1:ntotalbuckets
        push!(buckets, Vector{obj}())
    end
    return Buckets(buckets, worldbounds, nbuckets)
end

function isless_x(a::AABB, b::AABB)
    return (a.cornerpos[1] + a.sidelength[1]) < b.cornerpos[1]
end

minx(a::AABB) = a.cornerpos[1]
miny(a::AABB) = a.cornerpos[2]
maxx(a::AABB) = a.cornerpos[1] + a.sidelength[1]
maxy(a::AABB) = a.cornerpos[2] + a.sidelength[2]

#=
"""
Creates buckets in length(nbuckets) dimensions containing objects of type T
"""

function Buckets(obj::Type, worldaabb::AABB{2,T}, nbuckets::Tuple{Int,Int}) where {T}
    @assert length(nbuckets) == 2

    buckets = Vector{Bucket{2,T,obj}}()
    bucketpos = zeros(T, 2)
    
    bucketlength = (worldaabb.sidelength[1]./nbuckets[1], worldaabb.sidelength[2]./nbuckets[2])
    #This gives cryptic error
    #bucketlength = worldaabb.sidelength./nbuckets

    for ix in 1:nbuckets[1]
        bucketpos[1] = worldaabb.cornerpos[1] + (ix-1)*bucketlength[1]
        for iy in 1:nbuckets[2]
            bucketpos[2] = worldaabb.cornerpos[2] + (iy-1)*bucketlength[2]

            push!(buckets, Bucket(AABB(Tuple(bucketpos),Tuple(bucketlength)), obj[]))
        end
    end
    return Buckets(buckets, worldaabb, nbuckets)
end

function Buckets(obj::Type, worldaabb::AABB{3,T}, nbuckets::Tuple{Int,Int,Int}) where {T}
    @assert length(nbuckets) == 3

    buckets = Vector{Bucket{3,T,obj}}()
    bucketpos = zeros(T, 3)
    
    bucketlength = (worldaabb.sidelength[1]./nbuckets[1], worldaabb.sidelength[2]./nbuckets[2], worldaabb.sidelength[3]./nbuckets[3])
    #This gives cryptic error
    #bucketlength = worldaabb.sidelength./nbuckets
    for ix in 1:nbuckets[1]
        bucketpos[1] = worldaabb.cornerpos[1] + (ix-1)*bucketlength[1]
        for iy in 1:nbuckets[2]
            bucketpos[2] = worldaabb.cornerpos[2] + (iy-1)*bucketlength[2]
            for iz in 1:nbuckets[3]
                bucketpos[3] = worldaabb.cornerpos[3] + (iz-1)*bucketlength[3]

                push!(buckets, Bucket(AABB(Tuple(bucketpos),Tuple(bucketlength)), obj[]))
            end
        end
    end
    return Buckets(buckets, worldaabb, nbuckets)
end
=#

function isintersecting(X::T,Y::T,W::T,H::T, x::T, y::T, w::T, h::T) where T
    if((x >= X || x+w >= X) && (x <= X+W || x+w <= X+W))
        if((y >= Y || y+h >= Y) && (y <= Y+H || y+h <= Y+H))
            return true
        end
    end
    return false
end

function isintersecting(X::T,Y::T,Z::T,W::T,H::T,D::T, x::T, y::T, z::T, w::T, h::T, d::T) where T
    
    if((x >= X || x+w >= X) && (x <= X+W || x+w <= X+W))
        if((y >= Y || y+h >= Y) && (y <= Y+H || y+h <= Y+H))
            if((z >= Z || z+d >= Z) && (z <= Z+D || z+d <= Z+D))
                return true
            end
        end
    end
    return false

end

function isintersecting(rec1::AABB{dim,T}, rec2::AABB{dim,T}) where {dim,T}
    return isintersecting(rec1.cornerpos..., rec1.sidelength..., rec2.cornerpos..., rec2.sidelength...)
end

#=function insert!(buckets::Vector{Bucket{dim, T}}, object::T, pos) where {dim,T}

    for bucket in buckets
        if isintersecting(bucket.bounds..., pos..., po)
            push!(bucket.objects, object)
            return
        end
    end

    error("Warning: could not find bucket for object $object at pos ($(pos))")
end=#

function get_bucket_aabb(buckets::Buckets{3,T}, i::Int) where T
    xmin,ymin,zmin = buckets.worldbounds.cornerpos
    xmax,ymax,zmax = buckets.worldbounds.cornerpos + buckets.worldbounds.sidelength
    nx,ny,nz = buckets.nbuckets


end

function insert3!(buckets::Buckets{3,T}, object::obj, aabb::AABB{3,T}) where {T,obj}
    dim = 3
    corners = Vector{Vector{T}}()
    push!(corners, aabb.cornerpos + (aabb.sidelength.*(1,0,0)))
    push!(corners, aabb.cornerpos + (aabb.sidelength.*(0,0,1)))
    push!(corners, aabb.cornerpos + (aabb.sidelength.*(0,0,1)))
    push!(corners, aabb.cornerpos + (aabb.sidelength.*(1,1,0)))
    push!(corners, aabb.cornerpos + (aabb.sidelength.*(0,1,1)))
    push!(corners, aabb.cornerpos + (aabb.sidelength.*(1,0,1)))
    push!(corners, aabb.cornerpos + (aabb.sidelength.*(1,1,1)))
    push!(corners, aabb.cornerpos)


    xmin,ymin,zmin = buckets.worldbounds.cornerpos
    xmax,ymax,zmax = buckets.worldbounds.cornerpos + buckets.worldbounds.sidelength
    nx,ny,nz = buckets.nbuckets

    bb = Int[]
    sizehint!(bb,8)
    for ic in 1:8
        x,y,z = corners[ic]

        px = trunc(Int, nx*(x-xmin)/(xmax-xmin)+1) 
        py = trunc(Int, ny*(y-ymin)/(ymax-ymin)+1) 
        pz = trunc(Int, nz*(z-zmin)/(zmax-zmin)+1) 
        i = px + (py-1)*px + (pz-1)*px*py 

        if i > 0 || i < ((nx*ny*nz))
            if !(i in bb)
                push!(bb,i)
                push!(buckets.objects[i], object)
            end
        end
    end
end

function insert3!(buckets::Buckets{dim,T}, object::obj, pos::Vec{dim,T}) where {dim,T,obj}
    
    x,y,z = pos#80.0,10.0,30.0
    xmin,ymin,zmin = buckets.worldbounds.cornerpos
    xmax,ymax,zmax = buckets.worldbounds.cornerpos + buckets.worldbounds.sidelength
    nx,ny,nz = buckets.nbuckets

    px = trunc(Int, nx*(x-xmin)/(xmax-xmin)) + 1
    py = trunc(Int, ny*(y-ymin)/(ymax-ymin)) + 1
    pz = trunc(Int, nz*(z-zmin)/(zmax-zmin)) + 1
    
    i = px + (py-1)*px + (pz-1)*px*py

    if i > 0 && i <= (nx*ny*nz)
        push!(buckets.objects[i], object)
    end
end


function insert2!(buckets::Vector{Bucket{dim,T,obj}}, object::obj, aabb::AABB{dim,T}) where {dim,T,obj}
    
    n = false
    for bucket in buckets
        if isintersecting(bucket.bounds, aabb)
            n = true
            push!(bucket.objects, object)
            # return; dont call return since a mastersegment can fit into multiple buckets
        end
    end

    if n == false
        for bucket in buckets
            @show bucket
        end
        error("Warning: could not find bucket for object $object at pos $(aabb)")

    end
end


function empty_buckets!(buckets::Buckets)
	for bucket in buckets.buckets
		empty!(bucket.objects);
	end
end
