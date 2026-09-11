using BellMeasurementDynamics
using BenchmarkTools

const N = 100
const T = 200*N
const trials = 1_000
const p = 0.5
const alpha = 0.5


println("Threads: $(Threads.nthreads())")
println("N = $N")
println("T = $T")
println("trials = $trials")

@btime simulate_multiple($N, $T, $p, $alpha, $trials; threaded=true)