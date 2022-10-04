


module BeamExample
    mktempdir() do dir
        cd(dir) do
            include(joinpath(@__DIR__, "../docs/src/literate/beam_example.jl"))
        end
    end
end

module BarExample
    mktempdir() do dir
        cd(dir) do
            include(joinpath(@__DIR__, "../docs/src/literate/bar_example.jl"))
        end
    end
end

module EnfExample
    mktempdir() do dir
        cd(dir) do
            include(joinpath(@__DIR__, "../docs/src/literate/enf_example.jl"))
        end
    end
end

module PhaseFieldExample
    mktempdir() do dir
        cd(dir) do
            include(joinpath(@__DIR__, "../docs/src/literate/phasefield.jl"))
        end
    end
end