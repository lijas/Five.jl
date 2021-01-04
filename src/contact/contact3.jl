#https://www.osti.gov/servlets/purl/10175733

@with_kw mutable struct ContactSingleSurface{dim,T,CE<:AbstractContactEntity} <: AbstractContactSearchAlgorithm
	
	#slaves::Vector{AbstractContactEntity}
	#masters::Vector{FaceContactEntity{2,T,Lagrange{1,RefCube,1},2}}
    #masters::Vector{FaceContactEntity{3,T,Lagrange{2,RefCube,1},4}}
    masters::Vector{CE}
    bounding_boxes::Vector{AABB{dim,T}}

    neighbour_matrix::SparseMatrixCSC{Bool,Int64}

    nlbox::Vector{Int} #number of buckets a entity occupies
    lbox::Vector{Int} 
    nbox::Vector{Int} #number of entities in each bucket
    ndsort::Vector{Int} #list of entities sorted by bucket id
    npoints::Vector{Int} #pointer that identifies first entity in ndsort

    minx::T
    maxx::T

    miny::T
    maxy::T

    nbucketsx::Int
    nbucketsy::Int
    nbuckets::Int
    bucket_size::T

    possible_contacts::Vector{Tuple{Int,Int}}
end

function ContactSingleSurface{dim,T,CE}(master_segments::Vector{CE}) where {dim,T,CE<:AbstractContactEntity}
    nmasters = length(master_segments)
    bounding_boxes = Vector{AABB{dim,T}}(undef, nmasters)

    neighbour_matrix = spzeros(Bool,nmasters, nmasters)
    for i in 1:nmasters
        for j in i:nmasters

            share_node = false
            for ent_a_dof in master_segments[i].dofs
                for ent_b_dof in master_segments[j].dofs
                    if ent_a_dof[1] == ent_b_dof[1]

                        share_node = true
                    end
                end
            end
            if share_node == true
                neighbour_matrix[i,j] = true

            end
        end
    end

    return ContactSingleSurface{dim,T,CE}(master_segments, bounding_boxes, neighbour_matrix, zeros(Int, nmasters), 
                                        Int[], Int[], Int[], Int[], 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0.0,
                                        Tuple{Int,Int}[])
end

function update_contact!(contact::ContactSingleSurface{2,T}, x, it) where {T}

    if it%5 != 0
        return
    end
    dim = 2
    nmasters = length(contact.masters)

    @unpack_ContactSingleSurface contact

    minx = maxx = miny = maxy = 0.0
    bucket_size = Inf
    #Get maximum and minimum coordinates
    @timeit "Creating AABB" for (i, master) in enumerate(contact.masters)

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
    lbox = Int[]
    
    resize!(npoints, nbuckets)

    @timeit "Creating nlbox" for (i, master) in enumerate(contact.masters)

        aabb = bounding_boxes[i]

        min_coords = aabb.cornerpos
        max_coords = aabb.cornerpos + aabb.sidelength

        @timeit "trunc" begin
        xl = trunc(Int, (min_coords[1] - minx)/bucket_size) + 1
        xu = trunc(Int, (max_coords[1] - minx)/bucket_size) + 1
        
        yl = trunc(Int, (min_coords[2] - miny)/bucket_size) + 1
        yu = trunc(Int, (max_coords[2] - miny)/bucket_size) + 1
        end

        for ix in xl:xu
            for iy in yl:yu
                ib = (iy-1)*nbucketsx + ix
                nbox[ib] +=1
                push!(lbox,ib)
                nlbox[i] += 1
            end
        end
    end

    npoints[1] = 1
    @timeit "Npoints" for j in 2:nbuckets
        npoints[j] = npoints[j-1] + nbox[j-1]
    end
    resize!(ndsort, sum(nbox)+1)

    #@show nbox
    #@show lbox
    #@show nlbox
    #sort
    fill!(nbox, 0)
    c = 0
    @timeit "Sort contacts" for i in 1:nmasters
        for j in 1:nlbox[i]
            c += 1
            ib = lbox[c] #boxid
            ndsort[nbox[ib] + npoints[ib]] = i
            nbox[ib] += 1
        end
    end

    contact.maxx = maxx
    contact.maxy = maxy

    contact.minx = minx
    contact.miny = miny
    contact.bucket_size = bucket_size

    contact.nbucketsx = nbucketsx
    contact.nbucketsy = nbucketsy
    contact.nbuckets = nbuckets

    if maxx > 1000
        @show maxx, nbuckets
    end
    if minx < -1000
        @show minx, nbuckets
    end

    #ib = 1
    #@show nbox[ib]
    #@show ndsort[ (0:nbox[ib]) .+ npoints[ib]]
    #masterids in bucket ib
    #master_ids = ndsort[ (0:nbox[ib]) .+ npoint[ib] ]
    contact.possible_contacts = Tuple{Int,Int}[]
    @timeit "Possible contacts" for (i, master) in enumerate(masters)
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
                for j in 1:nbox[ib]
                    entity_id = ndsort[pointer + j-1]
                    push!(contact.possible_contacts, (i, entity_id))
                    #search_contact(master, masters[entity_id], x)
                end

            end
        end

    end

end

function update_contact!(contact::ContactSingleSurface{3,T}, x, it) where {T}

    if it%5 != 0
        return
    end
    dim = 3
    nmasters = length(contact.masters)

    @unpack_ContactSingleSurface contact

    #Instead of storing minx, miny etc... store them in min_dim = []
    minx = maxx = miny = maxy = 0.0
    min_dim = [Inf,Inf,Inf]
    max_dim = [-Inf,-Inf,-Inf]
    bucket_size = Inf
    #Get maximum and minimum coordinates
    @timeit "Creating AABB" for (i, master) in enumerate(contact.masters)

        aabb = getAABB(master, x)
        bounding_boxes[i] = aabb

        max_coords = aabb.cornerpos + aabb.sidelength
        for d in 1:dim
            min_dim[d] = min_dim[d] < aabb.cornerpos[d] ? min_dim[d] : aabb.cornerpos[d]
            max_dim[d] = max_dim[d] > max_coords[d] ? max_dim[d] : max_coords[d]
        end

        #bucketsize = is bases on the dimension of the smallest entity 
        max_dimension = ((aabb.sidelength[1]) > (aabb.sidelength[2])) ? (aabb.sidelength[1]) : (aabb.sidelength[2])
        if dim==3
            max_dimension = max_dimension > aabb.sidelength[3] ? max_dimension : aabb.sidelength[3]
        end
        bucket_size = (bucket_size > max_dimension) ? max_dimension : bucket_size
    end
    #bucket_size = 20.0
    nbuckets_dim = [0,0,0]
    nbuckets = 1
    for d in 1:dim
        nbuckets_dim[d] = trunc(Int, (max_dim[d]-min_dim[d])/bucket_size) + 1
        nbuckets *= nbuckets_dim[d]
    end

    resize!(nbox, nbuckets)
    fill!(nbox,0)
    fill!(nlbox,0)
    lbox = Int[]
    
    resize!(npoints, nbuckets)

    @timeit "Creating nlbox" for (i, master) in enumerate(contact.masters)

        aabb = bounding_boxes[i]

        min_coords = aabb.cornerpos
        max_coords = aabb.cornerpos + aabb.sidelength

        @timeit "trunc" begin
        xl = trunc(Int, (min_coords[1] - min_dim[1])/bucket_size) + 1
        xu = trunc(Int, (max_coords[1] - min_dim[1])/bucket_size) + 1
        
        yl = trunc(Int, (min_coords[2] - min_dim[2])/bucket_size) + 1
        yu = trunc(Int, (max_coords[2] - min_dim[2])/bucket_size) + 1
        end

        if dim == 2
            for ix in xl:xu
                for iy in yl:yu
                    ib = (iy-1)*nbucketsx + ix
                    nbox[ib] +=1
                    push!(lbox,ib)
                    nlbox[i] += 1
                end
            end
        elseif dim == 3
            zl = trunc(Int, (min_coords[3] - min_dim[3])/bucket_size) + 1
            zu = trunc(Int, (max_coords[3] - min_dim[3])/bucket_size) + 1
            for ix in xl:xu
                for iy in yl:yu
                    for iz in zl:zu
                        ib = (iz-1)*nbuckets_dim[1]*nbuckets_dim[2] + 
                             (iy-1)*nbuckets_dim[1] + ix

                        nbox[ib] +=1
                        push!(lbox,ib)
                        nlbox[i] += 1
                    end
                end
            end
        end
    end

    npoints[1] = 1
    @timeit "Npoints" for j in 2:nbuckets
        npoints[j] = npoints[j-1] + nbox[j-1]
    end
    resize!(ndsort, sum(nbox)+1)

    #@show nbox
    #@show lbox
    #@show nlbox
    #sort
    fill!(nbox, 0)
    c = 0
    @timeit "Sort contacts" for i in 1:nmasters
        for j in 1:nlbox[i]
            c += 1
            ib = lbox[c] #boxid
            ndsort[nbox[ib] + npoints[ib]] = i
            nbox[ib] += 1
        end
    end

    contact.maxx = max_dim[1]
    contact.maxy = max_dim[2]
    #contact.maxz = max_dim[3]

    contact.minx = min_dim[1]
    contact.miny = min_dim[2]
    #contact.minz = min_dim[3]
    contact.bucket_size = bucket_size

    contact.nbucketsx = nbuckets_dim[1]
    contact.nbucketsy = nbuckets_dim[2]
    #contact.nbucketsz = nbuckets_dim[3]
    contact.nbuckets = nbuckets
    #=
    @show nbuckets
    @show min_dim
    @show max_dim
    @show nbox
    @show bucket_size
    =#
    #ib = 1
    #@show nbox[ib]
    #@show ndsort[ (0:nbox[ib]) .+ npoints[ib]]
    #masterids in bucket ib
    #master_ids = ndsort[ (0:nbox[ib]) .+ npoint[ib] ]
    contact.possible_contacts = Tuple{Int,Int}[]
    contact_pairs = spzeros(Bool,nmasters, nmasters)
    @timeit "Possible contacts" for (i, master) in enumerate(masters)
        aabb = bounding_boxes[i]
        min_coords = aabb.cornerpos
        max_coords = aabb.cornerpos + aabb.sidelength

        ibox_min = min(nbuckets_dim[1], trunc(Int, (min_coords[1]-min_dim[1])/bucket_size)+1)
        ibox_max = min(nbuckets_dim[1], trunc(Int, (max_coords[1]-min_dim[1])/bucket_size)+1)

        jbox_min = min(nbuckets_dim[2], trunc(Int, (min_coords[2]-min_dim[2])/bucket_size)+1)
        jbox_max = min(nbuckets_dim[2], trunc(Int, (max_coords[2]-min_dim[2])/bucket_size)+1)

        if dim == 2
            for ix in ibox_min:ibox_max
                for iy in jbox_min:jbox_max
                    ib = (iy-1)*nbucketsx + ix

                    pointer = npoints[ib]
                    for j in 1:nbox[ib]
                        entity_id = ndsort[pointer + j-1]
                        push!(contact.possible_contacts, (i, entity_id))
                        #search_contact(master, masters[entity_id], x)
                    end

                end
            end
        elseif dim == 3
            kbox_min = min(nbuckets_dim[3], trunc(Int, (min_coords[3]-min_dim[3])/bucket_size)+1)
            kbox_max = min(nbuckets_dim[3], trunc(Int, (max_coords[3]-min_dim[3])/bucket_size)+1)

            for ix in ibox_min:ibox_max
                for iy in jbox_min:jbox_max
                    for iz in kbox_min:kbox_max
                        ib = (iz-1)*nbuckets_dim[1]*nbuckets_dim[2] + 
                             (iy-1)*nbuckets_dim[1] + ix

                        pointer = npoints[ib]

                        nbox[ib] == 0 ? continue : nothing #nothing to contact with
                        nbox[ib] == 1 ? continue : nothing #cant contact with itself

                        for j in 1:nbox[ib]
                            entity_id = ndsort[pointer + j-1]
                            low_id, high_id = minmax(i, entity_id)
                            if low_id == high_id
                                continue
                            end

                            is_neighbours = neighbour_matrix[low_id, high_id]
                            if is_neighbours == true
                                continue
                            end

                            if contact_pairs[low_id, high_id] == false
                                contact_pairs[low_id, high_id] = true
                                push!(contact.possible_contacts, (i, entity_id))
                            end
                        end
                    end

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

function search1!(contact::ContactSingleSurface, x, timestep, should_ouput = true)

    contact_outputs = []
    for contact_pair in contact.possible_contacts
        ent1 = contact.masters[contact_pair[1]]
        ent2 = contact.masters[contact_pair[2]]
        
        _co = search_contact(ent1, ent2, x)
        append!(_co, search_contact(ent2, ent1, x))
        for co in _co
            if co != nothing
                push!(contact_outputs, co)
            end
        end
    end
    return contact_outputs

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

function search_contact(ent_a::FaceContactEntity{dim,Float64,I,N}, ent_b::FaceContactEntity{dim,Float64,I,N}, x::Array{Float64,1}) where {dim,I,N}

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









