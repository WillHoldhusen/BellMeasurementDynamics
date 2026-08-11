@testset "alpha = 0" begin
    function test_unpaired(p)
        trials = 10_000
        N = 100
        alpha = 0.0
        T = 50*N

        n0 = 0.0
        S = 0.0

        for t in 1:trials
            state = simulate(N, T, p, alpha)
            n0 += BellMeasurementDynamics.unpaired_sites(state)/(N*trials)
            ents = BellMeasurementDynamics.entropy_profile(state)
            S += ents[div(N,2)]/trials 
        end
        a = sqrt(4*p-3*p^2)
        n0_expected = 0.5*(a-p)/(1-p)
        S_expected = 0.5*(1-p)*N/(2+a-p)


        @test n0 ≈ n0_expected atol=0.05

        @test S ≈ S_expected rtol=0.05

    end

    test_unpaired(0.1)
    test_unpaired(0.5)
end
