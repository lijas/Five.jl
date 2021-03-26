
# to test vtk-files
OVERWRITE_CHECKSUMS = false
checksums_file = joinpath(dirname(@__FILE__), "checksums.sha1")
checksum_list = read(checksums_file, String)
if OVERWRITE_CHECKSUMS
    csio = open(checksums_file, "w")
else
    csio = open(checksums_file, "r")
end

function get_MatLinearElastic_loading()
    mat = Five.MatLinearElastic(E = 200.0, nu = 0.3)
    
    strain1 = range(0.0, 1.0, length=100)
    strain2 = range(1.0, 0.5, length=100)
    strain3 = range(0.5, 2.0, length=100)

    _strain = [strain1..., strain2..., strain3...]
    strain = [SymmetricTensor{2,3}((x, 0.0, 0.0, 0.0, 0.0, 0.0)) for x in _strain]

    return (mat => strain)
end

function get_MatCZKolluri_loading()
    mat = MatCZKolluri(σₘₐₓ = 10.0, τₘₐₓ = 10.0, Φₙ = 1.0, Φₜ = 1.0)
    
    jump1 = collect(range(0.0,      stop = mat.δₙ*5, length=100))
    jump2 = collect(range(mat.δₙ*5, stop = 0.0, length=100))
    jump3 = collect(range(0.0,      stop = mat.δₙ*10, length=100))  

    _jump = [jump1..., jump2..., jump3...]
    jump = [Vec{3}((0.0, 0.0, x)) for x in _jump]

    return (mat => jump)
end

function get_MatCZKolluri_loading2()
    mat = MatCZKolluri(σₘₐₓ = 10.0, τₘₐₓ = 10.0, Φₙ = 1.0, Φₜ = 1.0)
    
    jump1 = collect(range(0.0,      stop = mat.δₜ*5, length=100))
    jump2 = collect(range(mat.δₜ*5, stop = 0.0, length=100))
    jump3 = collect(range(0.0,      stop = mat.δₜ*10, length=100))  

    _jump = [jump1..., jump2..., jump3...]
    jump = [Vec{3}((x/5, x/5, x/5)) for x in _jump]

    return (mat => jump)
end

function get_MatCZKolluri_loading3()
    mat = MatCZKolluri(σₘₐₓ = 10.0, τₘₐₓ = 10.0, Φₙ = 1.0, Φₜ = 1.0)
    
    jump1 = collect(range(0.0,      stop = mat.δₜ*5, length=100))
    jump2 = collect(range(mat.δₜ*5, stop = 0.0, length=100))
    jump3 = collect(range(0.0,      stop = mat.δₜ*10, length=100))  

    _jump = [jump1..., jump2..., jump3...]
    jump = [Vec{2}((x/2, x/2)) for x in _jump]

    return (mat => jump)
end

function get_MatCZKolluri_loading3()
    mat = MatCZKolluri(σₘₐₓ = 10.0, τₘₐₓ = 10.0, Φₙ = 1.0, Φₜ = 1.0)
    
    jump1 = collect(range(0.0,      stop = mat.δₜ*5, length=100))
    jump2 = collect(range(mat.δₜ*5, stop = 0.0, length=100))
    jump3 = collect(range(0.0,      stop = mat.δₜ*10, length=100))  

    _jump = [jump1..., jump2..., jump3...]
    jump = [Vec{2}((x/2, x/2)) for x in _jump]

    return (mat => jump)
end

@testset "Material behaviour" begin


    material_list = [
        get_MatLinearElastic_loading(),
        get_MatCZKolluri_loading(),
        get_MatCZKolluri_loading2(),
        get_MatCZKolluri_loading3()
    ]
        
    for (material, loading) in material_list
        state = Five.getmaterialstate(material)
        stresses = []; tangents = [];
        for load in loading
            stress, tangent, state = Five.constitutive_driver(material, load, state)
            push!(stresses, stress)
            push!(tangents, tangent)
        end

        checkhash1 = string(hash(stresses))
        checkhash2 = string(hash(tangents))
        if OVERWRITE_CHECKSUMS
            write(csio, checkhash1, "\n")
            write(csio, checkhash2, "\n")
        else
            @test chomp(readline(csio)) == checkhash1
            @test chomp(readline(csio)) == checkhash2
        end
    end

end

if OVERWRITE_CHECKSUMS
    close(csio)
end