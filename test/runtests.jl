using Test
using BellMeasurementDynamics

@testset "BellMeasurementDynamics" begin
    include("test_sampler.jl")
    include("test_alpha0.jl")
    include("test_threading.jl")
end