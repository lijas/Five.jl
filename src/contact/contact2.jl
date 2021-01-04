#https://www.osti.gov/servlets/purl/10175733

mutable struct ContactSingleSurface{dim,T} <: AbstractContactSearchAlgorithm
	
	#slaves::Vector{AbstractContactEntity}
	masters::Vector{AbstractContactEntity}

    slave_aabb::Vector{AABB{dim,T}}
    #master_aabb::Vector{AABB{dim,T}}

end

function ContactSingleSurface{dim,T}(master_segments::Vector{<:AbstractContactEntity}) where {dim,T}
    return ContactSingleSurface{dim,T}(master_segments, AABB{dim,T}[])
end

function update1!(contact::ContactSingleSurface{dim,T}, x, v) where {dim,T}

    nmasters = length(contact.masters)

    bounding_boxes = AABB{dim,T}[]
    Rx = zeros(Int, nmasters)
    Ry = zeros(Int, nmasters)

    for master in contact.masters
        aabb = getAABB(master, x)
        push!(bounding_boxes, aabb)
    end

    Ix = sortperm(bounding_boxes, by=(x)->x.cornerpos[1])
    Iy = sortperm(bounding_boxes, by=(x)->x.cornerpos[2])
    #Iz = sortperm(bounding_boxes, by=(x)->x.cornerpos[3])

    for i in 1:length(bounding_boxes)
        Rx[i] = Ix[i]
        Ry[i] = Iy[i]
        #Rz[i] = Iz[i]
    end

    for (i, master) in enumerate(contact.masters)
        range_x = searchsorted(Ix, bounding_boxes[i], by=(i)->bounding_boxes[Ix[i]], lt=isless_x)
        @show range_x

        _inside_box = Int[]
        for ix in range_x
            id = Ix[ix]
            if range_x[1] <= Ry[id] <= range_x[end]
                push!(_inside_box, id)
            end
        end

    end


end

function search1!(contact::ContactSingleSurface, x,v)

    for i = 1.2

    end

end