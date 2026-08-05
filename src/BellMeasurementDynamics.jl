module BellMeasurementDynamics
    export simulate, simulate_multiple, simulate_multiple_threaded
    include("AliasTables.jl")
    include("Simulation.jl")
end