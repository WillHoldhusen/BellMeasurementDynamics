@testset "threaded" begin
    N = 10
    T = 200*N
    p = 0.5
    alpha = 0.5
    trials = 100
    @test simulate_multiple(N, T, p, alpha, trials; threaded=false) ==
    simulate_multiple(N, T, p, alpha, trials; threaded=true)
end