export MatEGP

const PATH_TO_FORTRAN_SO = joinpath(@__DIR__,"egpmain")
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
    εᵖ::Float64
end

# # # # # # #
# Constructors
# # # # # # #

function getmaterialstate(material::MatEGP)
    
	dprops, iprops = getprops(material)

	ngens = iprops[8]
	nstats = iprops[7]
	_σ = zeros(Float64, ngens)
	_dσdε = zeros(Float64, ngens, ngens)	
	history = zeros(Float64, nstats)

    @assert(ngens == 6)
    @assert(length(dprops) == nstats)

    F0 = Matrix{Float64}(I, 3, 3)
	F1 = Matrix{Float64}(I, 3, 3)
    dprops[1] = 0.0
    
    _INIT_HISTORY_!(_σ,_dσdε, F0, F1, history, dprops, iprops)

    σ = fromvoigt(SymmetricTensor{2,3}, _σ, order = EO)

    return MatEGPState(σ, one(Tensor{2,3}), history, history[15*iprops[6]+2])
end

get_material_state_type(::MatEGP) = MatEGPState


# # # # # # #
# Drivers
# # # # # # #

function constitutive_driver(mp::MatEGP, F::Tensor{2,3}, state::MatEGPState, dt::Float64)
    ngens = mp.NGENS
    _σ = tovoigt(state.σ, order = EO)
    _dσdε = zeros(Float64, ngens, ngens)
    H = copy(state.H)

    #Set time step
    dprops, iprops = getprops(mp)
    dprops[1] = dt

    _CONST_DRIVER_!(_σ, _dσdε, Matrix(state.F), Matrix(F), H, dprops, iprops)

    σ = fromvoigt(SymmetricTensor{2,3}, _σ, order = EO)
    dσdε = fromvoigt(SymmetricTensor{4,3}, _dσdε, order = EO) 
    return σ, dσdε, MatEGPState(σ, F, H, H[15*iprops[6]+2])

    #J = det(F)
    #c = inv(J)*(otimesu(F,F) ⊡ dσdε ⊡ otimesu(F',F'))
    #S = J * inv(F) ⋅ σ ⋅ inv(F')
    #return S, c, MatEGPState(σ, F, H, H[15*iprops[6]+2])
end


function constitutive_driver(m::Material2D{MatEGP}, F::Tensor{2,2,T}, state::AbstractMaterialState, dt::Float64) where T
    @assert(m.plane_state == PLANE_STRAIN)

    #Convert to 3d
    F₃ = Tensor{2,3,T,9}((F[1,1], zero(T), F[1,2], zero(T), one(T), zero(T), F[2,1], zero(T), F[2,2]))
    σ₃, dσdε₃, newstate₃ = constitutive_driver(m.material, F₃, state, dt)

    #Convert back to 3d
    σ = SymmetricTensor{2,2,T,3}((σ₃[1,1],σ₃[1,3],σ₃[3,3]))
    dσdε = SymmetricTensor{4,2,T,9}((dσdε₃[1,1,1,1], dσdε₃[3,1,1,1], dσdε₃[3,3,1,1], dσdε₃[1,1,3,1], dσdε₃[3,1,3,1], dσdε₃[3,3,3,1], dσdε₃[1,1,3,3], dσdε₃[3,1,3,3], dσdε₃[3,3,3,3]))

    return σ, dσdε, newstate₃
end

# # # # # # #
# UTILS
# # # # # # #

function getprops(props::MatEGP)
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
    
    #dprops = zeros(Float64, iprops[9])
    dprops = [1.00000000e+01,2.73000000e+02,2.73000000e+02,2.60000000e+01,0.00000000e+00,0.00000000e+00,3.75000000e+03,1.00000000e-10,2.89000000e+05,0.00000000e+00,0.00000000e+00,0.00000000e+00,2.65000000e+01,9.65000000e-01,5.00000000e+01,-5.00000000e+00,5.38452923e-27,8.00000000e-02,3.21000000e+02,9.94474705e-45]
    #=dprops[1:8] .= props.STANDPROP

    start = zeros(Int, iprops[2])
    modes = zeros(Int, iprops[2])
    if iprops[2]==1
        start[1] = alphastart
        modes[1] = iprops[3]
    elseif iprops[2]==2
        start[1] = alphastart
        start[2] = alphastart+10+2+1*self.iprops[3]
        modes[1] = self.iprops[3]
        modes[2] = self.iprops[4]
    elseif iprops[2]==3
        start[1] = alphastart
        start[2] = alphastart+2*10+2+1*self.iprops[3]+2+1*self.iprops[4]
        start[3] = alphastart+2*10+2+1*self.iprops[3]+2+1*self.iprops[4]
        modes[1] = iprops[3]
        modes[2] = iprops[4]
        modes[3] = iprops[5]
    else
        error("NOP IS WRONGLY DEFINED")
    end
      
    kB = 1.38064852e-23
    @show kB
    R = 8.31e0

    @show props.PROPS_PROC1_STD
    k = 0
    for i in 1:iprops[1]
        @show start[i]
        dprops[start[i]:start[i]+modelength_std-1] = props.PROPS_PROC1_STD
        dprops[start[i]+modelength_std:start[i]+modelength_std+modes[i]-1] = props.PROPS_PROC1_G
        dprops[start[i]+modelength_std+modes[i]:start[i]+modelength_std+2*modes[i]-1] = props.PROPS_PROC1_GH0R
        
        if dprops[start[i]+9] > 1.0e-20
            dprops[start[i]+9] = kB * dprops[2]/dprops[start[i]+9] / 1.0e6
        end
            
        for l in 1:modes[i]
            @show start[i]+10+modes[i]+l
            self.dprops[start[i]+10+modes[i]+l] = dprops[start[i]+10+modes[i]+l] * exp(-1.0*dprops[start[i]]/R/dprops[2])
        end
    end=#

    return dprops, iprops
end

function _INIT_HISTORY_!(σ::Vector{Float64}, dσdε::Matrix{Float64}, F0::Matrix{Float64}, F1::Matrix{Float64}, history::Vector{Float64}, dprops::Vector{Float64}, iprops::Vector{Int32})
    iprops[1] = 6
	err = ccall((:egpgen_, PATH_TO_FORTRAN_SO), 
	Cvoid,
	(Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}), 
    σ, dσdε, F0, F1, history, dprops, iprops, 6|>Int32, length(history) |>Int32 , length(dprops) |> Int32)
end

function _CONST_DRIVER_!(σ::Vector{Float64}, dσdε::Matrix{Float64}, F0::Matrix{Float64}, F1::Matrix{Float64}, history::Vector{Float64}, dprops::Vector{Float64}, iprops::Vector{Int32})
    iprops[1] = 1
	err = ccall((:egpgen_, PATH_TO_FORTRAN_SO), 
	Cvoid,
	(Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}), 
	σ, dσdε, F0, F1, history, dprops, iprops, 6|>Int32, length(history) |>Int32 , length(dprops) |> Int32)
end