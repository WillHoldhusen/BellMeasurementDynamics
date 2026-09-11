module BellMeasurementDynamics
    export simulate!, simulate, simulate_multiple, average_entropy_profile, get_entropies_and_stats, simulate_multiple_time_series
    include("AliasTables.jl")
    include("Simulation.jl")
    include("Metrics.jl")
end