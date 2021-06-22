
#utility function for checking checksums of materials...
function check_checksum(material::Five.AbstractMaterial, loading::Vector, filename::String; OVERWRITE_CHECKSUMS = false)

    # to test vtk-files    
    checksums_file = joinpath(dirname(@__FILE__), "checksums", string(filename,".sha1"))
    #checksum_list = read(checksums_file, String)
    if OVERWRITE_CHECKSUMS
        csio = open(checksums_file, "w")
    else
        csio = open(checksums_file, "r")
    end

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

    if OVERWRITE_CHECKSUMS
        close(csio)
    end

end

include("test_MatCZKolluri.jl")
include("test_material2d.jl")
include("test_transv_mat.jl")
#include("test_MatHyperElasticPlastic.jl")