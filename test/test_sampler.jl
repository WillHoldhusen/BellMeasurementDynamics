@testset "Alias sampler" begin
    function test_sampler(weights)
        prob, alias = BellMeasurementDynamics.make_alias_table(weights)
        n = 1_000_000
        samples = BellMeasurementDynamics.sample(prob, alias, n)
        frequencies = [
            count(==(i), samples) / n
            for i in eachindex(weights)
        ]
        @test frequencies ≈ weights atol=0.005
    end

    test_sampler([.25, .25, .25, .25])
    test_sampler([.1, .2, .3, .4])
    test_sampler([.99, 0.0045, 0.005, 0.0005])

    alpha = 2.0
    N = 100
    r = 1:N
    real_weights = (N .- r).* (r.^(-alpha))
    real_weights = real_weights/sum(real_weights)
    test_sampler(real_weights)
end



