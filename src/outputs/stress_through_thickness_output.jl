export CellStressOutput 

struct CellStressOutput <: AbstractOutput
    cv::CellVectorValues
    center_quadpoints::Vector{Int}
end

function CellStressOutput(;part::AbstractPart{dim}) where dim

    @assert(part isa FEPart)
    @assert(part.element isa SolidElement)
    cv = copy(part.element.cv)

    #Get the quadrature points that are in the center of the element
    nqppoints = getnquadpoints(cv)
    _qporder = nqppoints^(1/dim)

    @assert( isinteger(_qporder) )
    qp_order = convert(Int, _qporder)
    @assert( isodd(qp_order) )

    row = ceil(Int, qp_order/2)
    startidx = qp_order*(row-1)*(dim-2) + row
    center_quadpoints = startidx:(qp_order^(dim-1)):nqppoints
    
    return CellStressOutput(cv, center_quadpoints)
end

function build_outputdata(output::CellStressOutput, set::Set{<:Int}, dh::MixedDofHandler)
    return output
end

function collect_output!(output::CellStressOutput, state::StateVariables, cellset::Set{Int}, globaldata)
    
    cv = output.cv

    tts = ThroughThicknessStresses()
    for cellid in sort(collect(cellset)) #NOTE, cellset should be a Vector to preserv order        
        Xe = cellcoords(globaldata.dh, cellid)

        reinit!(cv, Xe)

        materialstates = state.partstates[cellid].materialstates
        for qp in output.center_quadpoints

            R = _construct_rotaiton_matrix(cv, qp, Xe)
            _σ = materialstates[qp].σ 
            σ  = symmetric(R ⋅ _σ ⋅ R')

            push!(tts.stresses, σ)

            X = Ferrite.spatial_coordinate(cv, qp, Xe)

            xglob = Vec{3}((X[1],0.0,X[2]))
            push!(tts.global_coords, xglob)
            push!(tts.local_coords, R⋅xglob)
        end
    end

    return tts
end

function _construct_rotaiton_matrix(cv::CellVectorValues{dim,T}, qp::Int, X::Vector{Vec{dim,T}})::Tensor{2,3,T,9} where {dim,T}
    dxdξ = zeros(Vec{dim,T},dim-1)
    for i in 1:Ferrite.getngeobasefunctions(cv)
        dM = cv.dMdξ[i,qp]
        for d in 1:dim-1
            dxdξ[d] += dM[d] * X[i]
        end
    end

    dxdξ ./= norm.(dxdξ)
    D = Tensors.cross(dxdξ...)
    _R = hcat(dxdξ..., D)

    R = Tensor{2,dim}(Tuple(_R))

    if dim == 2
        return Tensor{2,3,T,9}((R[1,1], T(0.0), R[1,2], T(0.0), T(1.0), T(0.0), R[2,1], T(0.0), R[2,2]))
    end

    return R
end