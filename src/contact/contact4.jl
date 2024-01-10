#https://www.osti.gov/servlets/purl/10175733

struct NodeToSurfaceContact{dim,T,CE<:AbstractContactEntity} <: AbstractContactSearchAlgorithm
	
    masters::Vector{CE}
    nodes::Vector{NodeContactEntity{dim}}
    bounding_boxes::Vector{AABB{dim,T}}

    lbox::Vector{Int} 
    nbox::Vector{Int} #number of entities in each bucket
    ndsort::Vector{Int} #list of entities sorted by bucket id
    npoints::Vector{Int} #pointer that identifies first entity in ndsort

    nbuckets::Int
    bucket_size::T

    possible_contacts::Vector{Tuple{Int,Int}}
    xi_cache::Vector{Vec{2,T}}
    co_cache::Vector{ContactOutput{NodeContactEntity{dim},CE,dim,T}}
end

function NodeToSurfaceContact{dim,T,CE}(slave_nodes::Vector{NodeContactEntity}, master_surface::Vector{ContactElement}) where {dim,T,CE<:AbstractContactEntity}
    nmasters = length(master_segments)
    bounding_boxes = Vector{AABB{dim,T}}(undef, nmasters)
    
    nodes = NodeContactEntity{dim}[]
    _added_dofs = Int[]
    node_segements = Vector{Vector{Int}}()

    for imaster in 1:nmasters
        for j in 1:length(master_segments[imaster].dofs)
            master_dofs = master_segments[imaster].dofs[j]
            master_dof = master_dofs[1]
            idx = findfirst((x) -> x==master_dof, _added_dofs)
            if idx == nothing
                push!(nodes, NodeContactEntity(master_dofs))
                push!(_added_dofs, master_dof)
                push!(node_segements, Int[imaster])
            else
                push!(node_segements[idx], imaster)
            end
        end
    end

    return FeSurface{dim,T,CE}(master_segments, nodes, bounding_boxes, node_segements, zeros(Int, nmasters), 
                                        Int[], Int[], Int[], 0.0, 0.0,
                                        Tuple{Int,Int}[],
                                        Vec{2,T}[],
                                        Vector{ContactOutput{NodeContactEntity{dim},CE,dim,T}}(undef, 100))
end

function update_global_contact_sort!(contact::FeSurface{dim,T}, x) where {dim, T}

    nsegments = length(contact.masters)
    nnodes = length(contact.nodes)

    min_dim     = [Inf,Inf,Inf]
    max_dim     = [-Inf,-Inf,-Inf]
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
    
    nbuckets_dim = [0,0,0]
    nbuckets = 1
    for d in 1:dim
        nbuckets_dim[d] = trunc(Int, (max_dim[d]-min_dim[d])/bucket_size) + 1
        nbuckets *= nbuckets_dim[d]
    end

    resize!(nbox, nbuckets)
    fill!(nbox,0)
    lbox = zeros(Int, nnodes)
    
    resize!(npoints, nbuckets)

    for (i, node) in enumerate(nodes)

        node_coord = x[node.dofs]

        ix = trunc(Int, (node_coord[1] - min_dim[1])/bucket_size) + 1
        iy = trunc(Int, (node_coord[2] - min_dim[2])/bucket_size) + 1

        if dim == 2
            ib = (iy-1)*nbuckets_dim[1] + ix
            nbox[ib] +=1
            lbox[i] = ib

        elseif dim == 3
            iz = trunc(Int, (node_coord[3] - min_dim[3])/bucket_size) + 1

            ib = (iz-1)*nbuckets_dim[1]*nbuckets_dim[2] + 
                 (iy-1)*nbuckets_dim[1] + ix

            nbox[ib] +=1
            lbox[i] = ib

        end
    end

    npoints[1] = 1
    @timeit "Npoints" for j in 2:nbuckets
        npoints[j] = npoints[j-1] + nbox[j-1]
    end
    resize!(ndsort, sum(nbox)+1)

    #sort
    fill!(nbox, 0)
    @timeit "Sort contacts" for i in 1:nnodes
        ib = lbox[i] 
        ndsort[nbox[ib] + npoints[ib]] = i
        nbox[ib] += 1
    end

    #
    contact.possible_contacts = Tuple{Int,Int}[]
    contact.xi_cache = Vec{2,T}[]
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

                    nbox[ib] == 0 ? continue : nothing #nothing to contact with

                    for j in 1:nbox[ib]
                        node_id = ndsort[pointer + j-1]
                        
                        if i in node_segements[node_id]
                            continue
                        end

                        push!(potential_contact, (i, node_id))

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

                        for j in 1:nbox[ib]
                            node_id = ndsort[pointer + j-1]
                            
                            if i in node_segements[node_id]
                                 continue
                            end

                            push!(contact.possible_contacts, (i, node_id))
                            push!(contact.xi_cache, zero(Vec{2,T}))
                        end
                    end

                end
            end
        end

    end
end

function search1!(contact::FeSurface, x, timestep, should_ouput = true)

    #contact_outputs = []
    ncontacts = 0
    for (i, contact_pair) in enumerate(contact.possible_contacts)
        segment = contact.masters[contact_pair[1]]
        node = contact.nodes[contact_pair[2]]

        @timeit "local search" iscontact, co = search_contact(node, segment, x, contact.xi_cache[i])
        
        if iscontact == true
            ncontacts += 1
            contact.co_cache[ncontacts] = co
            contact.xi_cache[i] = co.xi

            if ncontacts > 100
                error("More than 100 contacts, cache overflow :(")
            end 
        end

    end
    return contact.co_cache, ncontacts

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









