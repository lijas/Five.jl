
#Function that calculates shell interpolation X₀ + 0.5t \zeta * D   
function _shell_interpolation(cv, x::AbstractVector{Vec{dim_s,T}}, ξ::Vec{dim_s,T2}) where {dim_s,T,T2}
    #dim_s = 2
    dim_p = dim_s-1
    t = cv.thickness
    msip = cv.mid_surface_interpolation

    nbasefunctions = JuAFEM.getnbasefunctions(msip)
    
    dX₀dξ = zeros(Vec{dim_s,T2}, dim_p)
    X₀ = zero(Vec{dim_s,T2})
    for i in 1:nbasefunctions
        _ξ = Vec((ξ[1:dim_p]...))
        dMdξ, M = gradient(ξ -> JuAFEM.value(msip, i, ξ), _ξ, :all)
        X₀ += M*x[i]
        for d in 1:dim_p
            dX₀dξ[d] += dMdξ[d]*x[i]
        end
    end

    _D = Tensors.cross(dX₀dξ...)
    D = _D/norm(_D)

    return X₀ + ξ[dim_s] * 0.5t * D

end

@testset "ShellCellValues" begin

    dim = 3
    T = Float64
    thickness = 1.0

    qr = QuadratureRule{dim,RefCube}(2)
    ip_geom = Lagrange{dim-1, RefCube, 2}()
    ip = Lagrange{dim, RefCube, 1}()

    scv = Five.ShellCellValues(thickness, qr, ip, ip_geom)

    qpoint = 1

    coords = Five.JuAFEM.reference_coordinates(ip_geom) .+ [rand(Vec{dim-1,T})*0.1 for _ in 1:Five.JuAFEM.getnbasefunctions(ip_geom)]
    coords = [Vec((x...,0.0)) for x in coords]
    J1 = Five._shell_interpolation(scv, qpoint, coords)
    J2 = gradient( ξ -> _shell_interpolation(scv, coords, ξ), scv.qr.points[qpoint])

    @test J1 ≈ J2
end