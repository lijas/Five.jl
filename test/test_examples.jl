


module BeamExample
    mktempdir() do dir
        cd(dir) do
            include(joinpath(@__DIR__, "../docs/src/literate/beam_example.jl"))
        end
    end
end

module BeamExample
    mktempdir() do dir
        cd(dir) do
            include(joinpath(@__DIR__, "../docs/src/literate/bar_example.jl"))
        end
    end
end