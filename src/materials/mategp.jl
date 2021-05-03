export MatEGP

const PATH_TO_FORTRAN_SO = "/home/elias/Documents/julia_code/Five.jl/src/materials/callFromJ"
const EO = [1 4 6; 7 2 5; 9 8 3] #Julia to EGP voigt order

"""
    MatEGP() <: AbstractMaterial

Eindoven glassy polymer material model
"""

Base.@kwdef struct MatEGP <: AbstractMaterial
    E::Float64 = 1.0e6
    ν::Float64 = 0.25
    
    ACTION::Int = 1;
    NOP::Int = 1
    NAM::Int = 1
    NBM::Int = 0
    NGM::Int = 0
    NGENS::Int = 6
    PRESMET::Int = 1
    HARDMET::Int = 1
    VISCHARD::Int = 0
    STIFFNESS::Int = 0
    VISCDEF::Int = 0

    STANDPROP::Vector{Float64} = [1.0E-2,2.73E2,2.73E2,2.6E1,0.0E0,0.0E0,3.75E3,1.0E0]

    PROPS_PROC1_STD::Vector{Float64} = [2.89E5,0.0E0,0.0E0,0.0E0,2.65E1,9.65E-1,5.0E1,-5.0E0,5.3845E-27,8.0E-2]
    PROPS_PROC1_G::Vector{Float64} = [3.21E2]
    PROPS_PROC1_GH0R::Vector{Float64} = [2.1E11]
end

struct MatEGPState <: AbstractMaterialState
    σ::SymmetricTensor{2,3,Float64,6}
    F::Tensor{2,3,Float64,9}
    H::Vector{Float64}
end

# # # # # # #
# Constructors
# # # # # # #

function getmaterialstate(material::MatEGP)
    
	iprops = getiprops(material)
	dprops = getdprops(material)

	ngens = iprops[8]
	nstats = iprops[7]
	_σ = zeros(Float64, ngens)
	_dσdε = zeros(Float64, ngens, ngens)	
	history = zeros(Float64, nstats)

    @assert(ngens == 6)
    @assert(length(dprops) == nstats)

    F0 = Matrix{Float64}(I, 3, 3)
	F1 = Matrix{Float64}(I, 3, 3)

    _INIT_HISTORY_!(_σ,_dσdε, F0, F1, history, dprops, iprops)

    σ = fromvoigt(SymmetricTensor{2,3}, _σ, order = EO)

    return MatEGPState(σ, one(Tensor{2,3}), history)
end

get_material_state_type(::MatEGP) = MatEGPState


# # # # # # #
# Drivers
# # # # # # #

function constitutive_driver(mp::MatEGP, F::Tensor{2,3}, state::MatEGPState)
    ngens = mp.NGENS
    _σ = zeros(Float64, ngens)
    _dσdε = zeros(Float64, ngens, ngens)
    H = copy(state.H)

    _CONST_DRIVER_!(_σ, _dσdε, Matrix(state.F), Matrix(F), H, getdprops(mp), getiprops(mp))

    σ = fromvoigt(SymmetricTensor{2,3}, _σ, order = EO)
    dσdε = fromvoigt(SymmetricTensor{4,3}, _dσdε, order = EO) 

    return σ, dσdε, MatEGPState(σ, F, H)
end

# # # # # # #
# UTILS
# # # # # # #

function getiprops(props::MatEGP)
	modelength_std = 10
    alphastart = 9

    # Read integer properies and store in self.iprops
    iprops = zeros(Int32, 14)
    iprops[1] = props.ACTION
    iprops[2] = props.NOP
    iprops[3] = props.NAM
    iprops[4] = props.NBM
    iprops[5] = props.NGM
    iprops[6] = props.NAM+props.NBM+props.NGM
    iprops[7] = 15*iprops[6] + 2*iprops[2] + 3
    iprops[8] = props.NGENS
    iprops[9] = alphastart-1 + iprops[2]*modelength_std + 2*iprops[6]
    iprops[10] = props.PRESMET
    iprops[11] = props.HARDMET
    iprops[12] = props.VISCHARD
    iprops[13] = props.STIFFNESS
    iprops[14] = props.VISCDEF
	return iprops
end

function getdprops(props::MatEGP)
	dprops = Float64[]
	append!(dprops, props.STANDPROP, props.PROPS_PROC1_STD, props.PROPS_PROC1_G, props.PROPS_PROC1_GH0R)
end

function _INIT_HISTORY_!(σ::Vector{Float64}, dσdε::Matrix{Float64}, F0::Matrix{Float64}, F1::Matrix{Float64}, history::Vector{Float64}, dprops::Vector{Float64}, iprops::Vector{Int32})
    iprops[1] = 6
	err = ccall((:egpgen_, PATH_TO_FORTRAN_SO), 
	Cvoid,
	(Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}), 
	σ, dσdε, F0, F1, history, dprops, iprops)
end

function _CONST_DRIVER_!(σ::Vector{Float64}, dσdε::Matrix{Float64}, F0::Matrix{Float64}, F1::Matrix{Float64}, history::Vector{Float64}, dprops::Vector{Float64}, iprops::Vector{Int32})
    iprops[1] = 1
	err = ccall((:egpgen_, PATH_TO_FORTRAN_SO), 
	Cvoid,
	(Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}), 
	σ, dσdε, F0, F1, history, dprops, iprops)
end