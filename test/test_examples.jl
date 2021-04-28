


module BeamExample
    mktempdir() do dir
        cd(dir) do
            include(joinpath(@__DIR__, "../docs/src/literate/beam_example.jl"))
        end
    end
end

#=module CurvedCantileverBeam
    using Test
    mktempdir() do dir
        cd(dir) do
            include(joinpath(@__DIR__, "../examples/stress_rec/curved_cantilever_beam/curved_cantilever_beam.jl"))
        end
    end
end=#