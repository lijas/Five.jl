#https://www.osti.gov/servlets/purl/10175733

mutable struct FeSurface{dim,T} <: AbstractContactSearchAlgorithm
	
	masterfaceset::Vector{FaceIndex}
    slavenodes   ::Vector{VertexIndex}

    masters::Vector{<:AbstractContactEntity}
    nodes::Vector{NodeContactEntity{dim}}

    bounding_boxes::Vector{AABB{dim,Float64}}

    node_segements::Vector{Vector{Int}} # node_segments[i] is a list of master elements containing node i

    nlbox::Vector{Int} #number of buckets a entity occupies
    lbox::Vector{Int}  #bucket id for each slave node
    nbox::Vector{Int} #number of entities in each bucket
    ndsort::Vector{Int} #list of entities sorted by bucket id
    npoints::Vector{Int} #pointer that identifies first entity in ndsort

    nbuckets::Int
    bucket_size::T

    possible_contacts::Vector{Tuple{Int,Int}}
end

function FeSurface{dim,T}(slavenodes, masterfaceset) where {dim,T}
    nmasters = length(masterfaceset)
    bounding_boxes = Vector{AABB{dim,T}}(undef, nmasters)
    
    nodes = NodeContactEntity{dim}[]
    master_segments = FaceContactEntity[]
    node_segements = Vector{Vector{Int}}()

    return FeSurface{dim,T}(masterfaceset, slavenodes,
                            master_segments, nodes, bounding_boxes, node_segements,
                            zeros(Int, nmasters), Int[], Int[], Int[], Int[], 
                            0, 0.0,
                            Tuple{Int,Int}[])
end

function Ferrite.close!(contact::FeSurface{dim,T}, dh::MixedDofHandler) where {dim,T}

    (; masterfaceset, slavenodes) = contact
    
    fh = getfieldhandler(dh, cellid(first(masterfaceset)))
    @assert(length(fh.fields)==1)
    faceinterp = Ferrite.getlowerdim(fh.fields[1].interpolation)
    N = getnbasefunctions(faceinterp)
    ip = typeof(faceinterp)

    masters = contact.masters #FaceContactEntity{dim,ip,N}[]
    slaves  = contact.nodes   #NodeContactEntity{dim}[]

    for faceidx in masterfaceset
        dofs = dofs_on_face(dh, fh, faceidx, :u, collect(1:dim))
        dofs = Tuple(reinterpret(NTuple{dim,Int},dofs))
        fce = FaceContactEntity{dim,ip,N}(faceinterp, dofs)
        push!(masters, fce)
    end

    for vertexidx in slavenodes
        dofs = dofs_on_vertex(dh, fh, vertexidx, :u, collect(1:dim))
        dofs = ntuple(i->dofs[i], dim)
        slavenode = NodeContactEntity{dim}(dofs)
        push!(slaves, slavenode)
    end
    
end

function update_contact!(contact::FeSurface{2,T}, x) where {T}

    (; masters, nodes, bounding_boxes, node_segements)     = contact
    (; lbox, nlbox, nbox, ndsort, npoints, nbuckets, bucket_size) = contact
    (; possible_contacts)              = contact

    dim = 2
    nmasters = length(contact.masters)

    minx = Inf
    maxx = -Inf
    miny = Inf
    maxy = -Inf
    bucket_size = -Inf
    #Get maximum and minimum coordinates
    for (i, master) in enumerate(contact.masters)

        aabb = getAABB(master, x)
        bounding_boxes[i] = aabb
        minx = minx < aabb.cornerpos[1] ? minx : aabb.cornerpos[1]
        miny = miny < aabb.cornerpos[2] ? miny : aabb.cornerpos[2]
        #minz = minz < aabb.cornerpos[3] ? minz : aabb.cornerpos[3]

        max_coords = aabb.cornerpos + aabb.sidelength
        maxx = maxx > max_coords[1] ? maxx : max_coords[1]
        maxy = maxy > max_coords[2] ? maxy : max_coords[2]
        #maxz = maxz > max_coords[3] ? maxz : max_coords[3]

        #bucketsize = is bases on the dimension of the smallest entity 
        max_dimension = ((aabb.sidelength[1]) > (aabb.sidelength[2])) ? (aabb.sidelength[1]) : (aabb.sidelength[2])
        bucket_size = (bucket_size > max_dimension) ? max_dimension : bucket_size
    end

    nbucketsx = trunc(Int, (maxx-minx)/bucket_size) + 1
    nbucketsy = trunc(Int, (maxy-miny)/bucket_size) + 1
    nbuckets = nbucketsy*nbucketsx

    resize!(nbox, nbuckets)
    fill!(nbox,0)
    fill!(nlbox,0)
    resize!(lbox, length(contact.nodes))
    
    resize!(npoints, nbuckets)

    for (i, node) in enumerate(contact.nodes)

        node_coord = x[Vec(node.dofs)]

        ix = trunc(Int, (node_coord[1] - minx)/bucket_size) + 1
        iy = trunc(Int, (node_coord[2] - miny)/bucket_size) + 1
        ib = (iy-1)*nbucketsx + ix
        nbox[ib] +=1
        lbox[i] = ib
    end

    npoints[1] = 1
    for j in 2:nbuckets
        npoints[j] = npoints[j-1] + nbox[j-1]
    end
    resize!(ndsort, sum(nbox)+1)

    fill!(nbox, 0)
    for i in 1:length(contact.slavenodes)
        ib = lbox[i] #boxid
        ndsort[nbox[ib] + npoints[ib]] = i
        nbox[ib] += 1
    end

    contact.nbuckets = nbuckets

    if maxx > 1000
        @show maxx, nbuckets
    end
    if minx < -1000
        @show minx, nbuckets
    end

    #ib = 1
    #@show ndsort[ (0:nbox[ib]) .+ npoints[ib]]
    #masterids in bucket ib
    #master_ids = ndsort[ (0:nbox[ib]) .+ npoint[ib] ]
    empty!(contact.possible_contacts)
    for (i, master) in enumerate(masters)
        aabb = bounding_boxes[i]
        min_coords = aabb.cornerpos
        max_coords = aabb.cornerpos + aabb.sidelength

        ibox_min = min(nbucketsx, trunc(Int, (min_coords[1]-minx)/bucket_size)+1)
        ibox_max = min(nbucketsx, trunc(Int, (max_coords[1]-minx)/bucket_size)+1)

        jbox_min = min(nbucketsy, trunc(Int, (min_coords[2]-miny)/bucket_size)+1)
        jbox_max = min(nbucketsy, trunc(Int, (max_coords[2]-miny)/bucket_size)+1)

        for ix in ibox_min:ibox_max
            for iy in jbox_min:jbox_max
                ib = (iy-1)*nbucketsx + ix

                pointer = npoints[ib]
                nbox[ib] == 0 ? continue : nothing #nothing to contact with

                for j in 1:nbox[ib]
                    node_id = ndsort[pointer + j-1]
                    #if i in node_segements[node_id]
                    #    continue
                    #end
                    push!(contact.possible_contacts, (i, node_id))
                end

            end
        end

    end
end

#=struct FiniteElementMesh{dim,T} <:AbstractContactEntity
    nodes::Vector{Vec{dim,T}}
    elements::Vector{Cells}

    bucket_sort:
end=#

function search1!(contact::FeSurface, x)
    #contact_outputs = []
    ncontacts = 0
    contacts = Any[]
    for (i, contact_pair) in enumerate(contact.possible_contacts)
        segment = contact.masters[contact_pair[1]]
        node = contact.nodes[contact_pair[2]]

        @timeit "local search" iscontact, co = search_contact(node, segment, x)#, contact.xi_cache[i])
        if iscontact == true
            ncontacts += 1
            push!(contacts, co)
            #contact.co_cache[ncontacts] = co
            #contact.xi_cache[i] = co.xi
            if ncontacts > 100
                error("More than 100 contacts, cache overflow :(")
            end 
        end

    end
    return contacts, ncontacts

end

#=function search_contact(ent_a::FaceContactEntity{2,Float64,Lagrange{1,RefCube,1},2}, ent_b::FaceContactEntity{2,Float64,Lagrange{1,RefCube,1},2}, x::Array{Float64,1})

    x1 = x[ent_a.dofs[1]]
    x2 = x[ent_a.dofs[2]]

    contact_outputs = []
    for i in 1:2
        #xs = x[ent_b.dofs[i]]
        if ent_b.dofs[i] == ent_a.dofs[1] || ent_b.dofs[i] == ent_a.dofs[2]
            continue
        end
        ne = NodeContactEntity(ent_b.dofs[i])
        co =search_contact(ne, ent_a, x, false)

        #co = _search_node_line_contact(xs,x1,x2)
        if co != nothing
            push!(contact_outputs, co)
        end
    end
    return contact_outputs

end=#

function search_contact(ent_a::FaceContactEntity{dim,I,N}, ent_b::FaceContactEntity{dim,I,N}, x::Array{Float64,1}) where {dim,I,N}

    contact_outputs = []
    for i in 1:N

        ne = NodeContactEntity(ent_b.dofs[i])
        co =   search_contact(ne, ent_a, x, false)

        if co != nothing
            push!(contact_outputs, co)
        end
    end
    return contact_outputs

end









