module FiveImmersed

using Five
using FerriteImmersed: GhostInterface, ImmersoCellValues
using IGA

struct ImmersoPartCache{dim,CV,GIV,M}
    ue::Vector{Float64}
    due::Vector{Float64}
    Δue::Vector{Float64}
    fe::Vector{Float64}
    ke::Matrix{Float64}
    celldofs::Vector{Int}
    bcoords::Vector{Vec{dim,Float64}}
    coords::Vector{Vec{dim,Float64}}
    bw::Vector{Float64}
    w::Vector{Float64}
    cv::CV
    δɛ::Vector{SymmetricTensor{2,dim,Float64,M}}
    materialoptions::Dict{Symbol, Any} #Yuck

    #For ghostinterface
    giv::GIV
    fe_penalty::Vector{Float64}
    ke_penalty::Matrix{Float64}
    ue_penalty::Vector{Float64}
    interfacedofs::Vector{Int}
    celldofs2::Vector{Int}
    bcoords2::Vector{Vec{dim,Float64}}
    coords2::Vector{Vec{dim,Float64}}
    bw2::Vector{Float64}
    w2::Vector{Float64}
end

function ImmersoPartCache{dim}(ndofs, ncoords, cv_template::Ferrite.Values, giv_template::Ferrite.Values) where dim
    T = Float64
    ImmersoPartCache(
        zeros(T,ndofs), 
        zeros(T,ndofs), 
        zeros(T,ndofs), 
        zeros(T,ndofs), 
        zeros(T,ndofs,ndofs), 
        zeros(Int,ndofs), 
        zeros(Vec{dim,T},ncoords),
        zeros(Vec{dim,T},ncoords),
        zeros(T,ncoords), 
        zeros(T,ncoords), 
        deepcopy(cv_template),
        zeros(SymmetricTensor{2,dim,Float64}, ndofs),
        Dict{Symbol, Any}(),
        
        deepcopy(giv_template),
        zeros(T,ndofs*2), 
        zeros(T,ndofs*2,ndofs*2), 
        zeros(T,ndofs*2), 
        zeros(Int,ndofs*2), 
        zeros(Int,ndofs), 
        zeros(Vec{dim,T},ncoords),
        zeros(Vec{dim,T},ncoords),
        zeros(T,ncoords), 
        zeros(T,ncoords)
    )
end

struct GhostPenalty
    γ::Float64
    interfaces::Vector{GhostInterface}
end

#Name it Part instead of Fe-part because it is the standard...
struct ImmersoPart{dim,MAT<:Five.MaterialModels.AbstractMaterial,IP<:IGA.Bernstein,CACHE<:ImmersoPartCache} <: Five.AbstractPart{dim}
    material::MAT
    cellset::Vector{Int}
    threadsets::Vector{Vector{Int}}
    ip::IP
    domainpos::NTuple{dim,Float64}
    domainsize::NTuple{dim,Float64}
    nvoxels_dim::NTuple{dim,Int}
    quadraturerules::Dict{Int,QuadratureRule{dim,RefCube,Float64}}
    ghost_penalty::Five.Optional{GhostPenalty}
    cellid_ghostinterface_map::Matrix{GhostInterface}
    cache::Vector{CACHE} 
    geometry::Five.Optional{ImmersoGeometry{dim}} 
end

function ImmersoPart(;
    material::Five.AbstractMaterial,
    cellset,
    domainpos::NTuple{dim,Float64},
    domainsize::NTuple{dim,Float64},
    nvoxels_dim::NTuple{dim,Int},
    quadraturerules::Dict{Int,QuadratureRule{dim,RefCube,Float64}},
    ghost_penalty::Five.Optional{GhostPenalty} = nothing,
    geometry::Five.Optional{ImmersoGeometry} = nothing
    ) where {dim}

    T = Float64

    ip_geo   = Bernstein{dim,2}()
    field_ip = ip_geo^dim
    _set = sort!(collect(cellset))

    #Create cellvalues
    cv = FerriteImmersed.ImmersoCellValues(field_ip)
    cellsize = domainsize ./ nvoxels_dim
    reinit!(cv, cellsize)

    qr = InterfaceQuadratureRule{dim,RefCube}(3)
    gv = GhostInterfaceValues(qr, field_ip)

    ncelldofs = getnbasefunctions(field_ip)
    ncellnodes = getnbasefunctions(ip_geo)
    
    nchunks = (Threads.nthreads()-1) * 10 + 1
    
    cache1 = ImmersoPartCache{dim}(ncelldofs, ncellnodes, cv, gv)
    cache  = [cache1]
    for i in 2:nchunks
        cache[i] = ImmersoPartCache{dim}(ncelldofs, ncellnodes, cv, gv)
    end

    local cellid_ghostinterface_map
    if ghost_penalty !== nothing
        #Check that the interfaces are part of the cellset
        @assert all(map(id -> (id.faceindex1[1] in _set), ghost_penalty.interfaces))
        @assert all(map(id -> (id.faceindex2[1] in _set), ghost_penalty.interfaces))
        #@assert mapreduce(&, id -> id.faceindex2[2] in cellset, ghost_penalty.interfaces, true)
        cellid_ghostinterface_map = Matrix{GhostInterface}(:undef, nfaces, length(_set))
        for interface in ghost_penalty.interfaces
            (cellid, faceid) = interface.faceindex1
            i = searchsortedfirst(_set, cellid)
            cellid_ghostinterface_map[faceid,i] = interface
        end
    else
        cellid_ghostinterface_map = Matrix{GhostInterface}(undef, 0, 0)
    end

    return ImmersoPart(
        material, 
        _set,
        Vector{Int}[],
        field_ip,
        domainpos,
        domainsize,
        nvoxels_dim,
        quadraturerules,
        ghost_penalty,
        cellid_ghostinterface_map,
        cache,
        geometry)
end

struct ImmersoPartState <: Five.AbstractPartState

end

Five.get_fields(part::ImmersoPart{dim}) where dim = [(:u, part.ip),]
Five.get_cellset(part::ImmersoPart) = part.cellset

function Five.construct_partstates(part::ImmersoPart) 
    return ImmersoPartState()
end

function Five.init_part!(part::ImmersoPart{dim, T}, dh::Ferrite.AbstractDofHandler) where {dim,T}
    grid = dh.grid
    ip = part.ip
    cellset = part.cellset

    #Create cellvalues
    cv = ImmersNurbs.ImmersoCellValues(ip)
    cellsize = part.domainsize ./ part.nvoxels_dim
    reinit!(cv, cellsize)

    qr = InterfaceQuadratureRule{dim,RefCube}(3)
    gv = GhostInterfaceValues(qr, ip)

    ncelldofs = getnbasefunctions(cv)
    ncellnodes = getnbasefunctions(ip)
    
    nchunks = (Threads.nthreads()-1) * 10 + 1
    
    resize!(part.cache, nchunks)
    for i in 1:nchunks
        part.cache[i] = ImmersoPartCache{dim}(ncelldofs, ncellnodes, cv, gv)
    end

    threadsets = Ferrite.create_coloring(grid, part.cellset; alg=ColoringAlgorithm.WorkStream) 
    copy!(part.threadsets, threadsets)

    @assert dh.grid.cells[first(part.cellset)] isa IGA.BezierCell{dim,2}
    @assert issorted(cellset)
end

function Five.assemble_stiffnessmatrix_and_forcevector!(
    dh        ::Ferrite.AbstractDofHandler, 
    part      ::ImmersoPart{sdim},
    partstate ::ImmersoPartState,
    state     ::Five.StateVariables) where sdim

    nfaces = 2*sdim
    nchunks = length(part.cache)
    assemblers = [start_assemble(state.system_arrays.Kⁱ, state.system_arrays.fⁱ, fillzero=false) for _ in 1:nchunks]

    #Loop over cells
    for color in part.threadsets
        Threads.@threads for (chunkrange, ichunk) in chunks(color, nchunks) 
            for i in chunkrange
                cellid = color[i]
                cache = part.cache[ichunk]
                assembler = assemblers[ichunk]
                
                #1. Assemble stiffness matrix
                _assemble_cell(dh, part.quadraturerules, cellid, cache, part.material, assembler, state.d)

                #2. Ghost interface
                if part.ghost_penalty !== nothing
                    for lfaceid in 1:size(part.cellid_ghostinterface_map, 1)
                        j = searchsortedfirst(part.cellset, cellid)
                        interface = part.cellid_ghostinterface_map[lfaceid, j]
                        if interface.faceindex1[1] != -1
                            _assemble_ghost_penalty(dh, part, state, assembler, interface, part.ghost_penalty.γ, cache)
                        end
                    end
                end
            end
        end
    end
end

function _assemble_cell(dh::DofHandler{sdim}, quadraturerules, cellid, cache, material, assembler, d) where sdim
    (; fe, ke, celldofs, cv, ue, δɛ, materialoptions) = cache
    stresstate = sdim == 3 ? MaterialModels.Dim{3}() : MaterialModels.PlaneStrain()

    qr = quadraturerules[cellid]
    beo = get_extraction_operator(dh.grid, cellid)

    fill!(fe, 0.0)
    fill!(ke, 0.0)
    celldofs!(celldofs, dh, cellid)

    Five.disassemble!(ue, d, celldofs)

    #Build cellvalues for current cell
    ImmersNurbs.set_quadraturerule!(cv, qr)
    ImmersNurbs.set_bezier_operator!(cv, beo)

    _integrate_element!(ke, fe, stresstate, material, cv, δɛ, materialoptions, ue)

    #Assemble in to stiffness matrix
    assemble!(assembler, celldofs, ke, fe)
end

function _integrate_element!(ke::AbstractMatrix, fe::AbstractVector, stresstate, material::MaterialModels.AbstractMaterial, cv::ImmersNurbs.ImmersoCellVectorValues{sdim}, δɛ, materialoptions, ae::Vector{Float64}) where sdim
    n_basefuncs = getnbasefunctions(cv)
    
    zerostate =  MaterialModels.initial_material_state(material)
    cache = nothing
    options = materialoptions

    Ω = 0.0
    for q_point in 1:getnquadpoints(cv)
        ImmersNurbs.reinit_qp!(cv, q_point)
        
        ε = symmetric(function_gradient(cv, q_point, ae))
        σ, C, _ = MaterialModels.material_response(stresstate, material, ε, zerostate, nothing; cache, options)

        for i in 1:n_basefuncs
            δɛ[i] = symmetric(shape_gradient(cv, q_point, i)) 
        end

        dΩ = getdetJdV(cv, q_point)
        for col in 1:n_basefuncs
            fe[col] += (σ ⊡ δɛ[col]) * dΩ
            Cδɛ = C ⊡ δɛ[col]
            for row in col:n_basefuncs
                ke[row,col] += (δɛ[row] ⊡ Cδɛ) * dΩ
            end
        end
    end

    symmetrize_upper!(ke)
end;

function symmetrize_upper!(K::Matrix)
    for i in 1:size(K,1)
        for j in i+1:size(K,1)
            K[i,j] = K[j,i]
        end
    end
end;

function Five.post_part!(dh, part::ImmersoPart, states::Five.StateVariables)
    
end

function Five.commit_part!(dh, part::ImmersoPart, state::Five.StateVariables)
    return nothing
end

function _assemble_ghost_penalty(dh, part, state, assembler, interface::GhostInterface, γ::Float64, cache)

    (; celldofs, celldofs2, coords, coords2, bcoords, bcoords2, w, w2, bw, bw2,interfacedofs) = cache
    (; ke_penalty, fe_penalty, ue_penalty, giv) = cache

    n = length(celldofs)
    r1 = (1:n) .+ 0
    r2 = (1:n) .+ n

    k = 2
    h = first(part.domainsize./part.nvoxels_dim)
    _factor = γ * h^(2k+1)/h^2

    (cellid1, lfaceid1) = interface.faceindex1
    (cellid2, lfaceid2) = interface.faceindex2
    fill!(ke_penalty, 0.0)
    fill!(fe_penalty, 0.0)

    get_bezier_coordinates!(bcoords, bw, coords, w, dh.grid, cellid1)
    get_bezier_coordinates!(bcoords2, bw2, coords2, w2, dh.grid, cellid2)
    celldofs!(celldofs,  dh, cellid1)
    celldofs!(celldofs2, dh, cellid2)
    interfacedofs[r1] .= celldofs
    interfacedofs[r2] .= celldofs2
    #totalcelldofs = vcat(celldofs,celldofs2)
    #ue_penalty = state.d[totalcelldofs]

    Five.disassemble!(ue_penalty, state.d, interfacedofs)

    set_bezier_operator!(
        giv, 
        get_extraction_operator(dh.grid, cellid1), 
        get_extraction_operator(dh.grid, cellid2)
    )
    reinit!(giv, bcoords, lfaceid1, bcoords2, lfaceid2)

    L = 0.0
    for iqp in 1:getnquadpoints(giv)
        dV = getdetJdV(giv, iqp)
        L += dV
        ∇∇u_jump = function_jump_normal_second_derivative(giv, iqp, ue_penalty)
        for i in 1:getnbasefunctions(giv)
            Ni = ImmersNurbs.jump_in_shape_normal_second_derivative(giv, iqp, i)
            fe_penalty[i] += (∇∇u_jump ⋅ Ni) * dV
            for j in 1:getnbasefunctions(giv)
                Nj = ImmersNurbs.jump_in_shape_normal_second_derivative(giv, iqp, j)
                ke_penalty[i,j] += (Ni ⋅ Nj) * dV
            end
        end
    end
    ke_penalty .*= _factor
    fe_penalty .*= _factor
    
    _hot_fix_nonunique_dofs_assemble!(assembler, interfacedofs, ke_penalty, fe_penalty, false)
end

function _assemble_ghost_penalty(dh, part, state, assembler, ghost_interfaces::Vector, γ::Float64)

    (; celldofs, celldofs2, coords, coords2, bcoords, bcoords2, w, w2, bw, bw2,) = part.cache[1]
    (; ke_penalty, fe_penalty, giv) = part.cache[1]

   # assembler = start_assemble()

    k = 2
    h = first(part.domainsize./part.nvoxels_dim)
    _factor = γ * h^(2k+1)/h^2

    for interface in ghost_interfaces
        (cellid1, lfaceid1) = interface.faceindex1
        (cellid2, lfaceid2) = interface.faceindex2
        fill!(ke_penalty, 0.0)
        fill!(fe_penalty, 0.0)

        get_bezier_coordinates!(bcoords, bw, coords, w, dh.grid, cellid1)
        get_bezier_coordinates!(bcoords2, bw2, coords2, w2, dh.grid, cellid2)
        celldofs!(celldofs,  dh, cellid1)
        celldofs!(celldofs2, dh, cellid2)
        totalcelldofs = vcat(celldofs,celldofs2)
        ae = state.d[totalcelldofs]

        set_bezier_operator!(
            giv, 
            get_extraction_operator(dh.grid, cellid1), 
            get_extraction_operator(dh.grid, cellid2)
        )
        reinit!(giv, bcoords, lfaceid1, bcoords2, lfaceid2)

        L = 0.0
        for iqp in 1:getnquadpoints(giv)
            dV = getdetJdV(giv, iqp)
            L += dV
            ∇∇u_jump = function_jump_normal_second_derivative(giv, iqp, ae)
            for i in 1:getnbasefunctions(giv)
                Ni = ImmersNurbs.jump_in_shape_normal_second_derivative(giv, iqp, i)
                fe_penalty[i] += (∇∇u_jump ⋅ Ni) * dV
                for j in 1:getnbasefunctions(giv)
                    Nj = ImmersNurbs.jump_in_shape_normal_second_derivative(giv, iqp, j)
                    ke_penalty[i,j] += (Ni ⋅ Nj) * dV
                end
            end
        end
        ke_penalty .*= _factor
        fe_penalty .*= _factor
        
        #minike, minife, minidofs = mini_condense(ke_penalty, fe_penalty, totalcelldofs)
        #_assemble!(assembler, minidofs, minike, minife)
        #@show cellid1 cellid2 totalcelldofs
        _hot_fix_nonunique_dofs_assemble!(assembler, totalcelldofs, ke_penalty, fe_penalty, false)
        #assemble!(assembler, totalcelldofs, ke_penalty)
        #assemble!(assembler, sorteddofs, miniKe)
    end

   # a = assembler
   # K2 = Ferrite.SparseArrays.sparse(a.I, a.J, a.V, ndofs(dh), ndofs(dh))
   # @show norm(K2)
    #return K2
    #state.system_arrays.Kⁱ.nzval .+= 100*eps(Float64)
    #state.system_arrays.Kⁱ .+= K2

end


Base.@propagate_inbounds function _hot_fix_nonunique_dofs_assemble!(A::Ferrite.AbstractSparseAssembler, dofs::AbstractVector{Int}, Ke::AbstractMatrix, fe::AbstractVector, sym::Bool)
    ld = length(dofs)
    @assert size(Ke, 1) == ld
    @assert size(Ke, 2) == ld
    if length(fe) != 0
        @assert length(fe) == ld
        @boundscheck checkbounds(A.f, dofs)
        @inbounds assemble!(A.f, dofs, fe)
    end

    K = Ferrite.matrix_handle(A)
    permutation = A.permutation
    sorteddofs = A.sorteddofs
    @boundscheck checkbounds(K, dofs, dofs)
    resize!(permutation, ld)
    resize!(sorteddofs, ld)
    copyto!(sorteddofs, dofs)
    Ferrite.sortperm2!(sorteddofs, permutation)

    current_col = 1
    @inbounds for Kcol in sorteddofs
        maxlookups = sym ? current_col : ld
        Kecol = permutation[current_col]
        ri = 1 # row index pointer for the local matrix
        Ri = 1 # row index pointer for the global matrix
        nzr = nzrange(K, Kcol)
        while Ri <= length(nzr) && ri <= maxlookups
            R = nzr[Ri]
            Krow = K.rowval[R]
            Kerow = permutation[ri]
            val = Ke[Kerow, Kecol]
            if Krow == dofs[Kerow]
                # Match: add the value (if non-zero) and advance the pointers
                if !iszero(val)
                    K.nzval[R] += val
                end
                ri += 1
                #Ri += 1 FIVE REMOVED THIS LINE
            elseif Krow < dofs[Kerow]
                # No match yet: advance the global matrix row pointer
                Ri += 1
            else # Krow > dofs[Kerow]
                # No match: no entry exist in the global matrix for this row. This is
                # allowed as long as the value which would have been inserted is zero.
                iszero(val) || error("some row indices were not found")
                # Advance the local matrix row pointer
                ri += 1
            end
        end
        # Make sure that remaining entries in this column of the local matrix are all zero
        for i in ri:maxlookups
            if !iszero(Ke[permutation[i], Kecol])
                error("some row indices were not found")
            end
        end
        current_col += 1
    end
end


function Five.assemble_sparsity_pattern!(part::ImmersoPart, globaldata)
    if part.ghost_penalty === nothing
        return sparse(Int[], Int[], Float64[], ndofs(globaldata.dh), ndofs(globaldata.dh))
    end

    celldofs1 = part.cache[1].celldofs
    celldofs2 = part.cache[1].celldofs2
    I = Int[]
    J = Int[]
    for interface in part.ghost_penalty.interfaces
        (cellid1, lfaceid1) = interface.faceindex1
        (cellid2, lfaceid2) = interface.faceindex2
        celldofs!(celldofs1,  globaldata.dh, cellid1)
        celldofs!(celldofs2,  globaldata.dh, cellid2)
        global_dofs = vcat(celldofs1, celldofs2)
        for i in eachindex(global_dofs)
            for j in eachindex(global_dofs)
                dofi = global_dofs[i]
                dofj = global_dofs[j]
                push!(I, dofi)
                push!(J, dofj)
            end
        end
    end
    _ndofs = ndofs(globaldata.dh)
    K = Ferrite.spzeros!!(Float64, I, J, _ndofs, _ndofs)
    fill!(K.nzval, 1)
    return K
end


end #module