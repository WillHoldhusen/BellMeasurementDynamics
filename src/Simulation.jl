include("Metrics.jl")
include("AliasTables.jl")

function measure_one!(s::Int, state::Vector{Int})
    j = state[s]
    if j != 0
        state[s] = 0
        state[j] = 0
    end
end

function measure_two!(s1::Int, s2::Int, state::Vector{Int})
    j1 = state[s1]
    j2 = state[s2]
    if j1 == 0 && j2 == 0 # both sites unentangled
        # perform Bell measurement entangling sites s1, s2
        state[s1] = s2
        state[s2] = s1
    elseif j1 != 0 && j2 != 0 # both sites already entangled
        if j1 != s2 && j2 != s1 # not entangled together
            state[s1] = s2 # entangle s1 and s2
            state[s2] = s1
            state[j1] = j2 # entangle the previous partners
            state[j2] = j1
        end
    elseif j1 == 0 && j2 != 0
        state[s1] = s2
        state[s2] = s1
        state[j2] = 0
    elseif j1 != 0 && j2 == 0
        state[s1] = s2
        state[s2] = s1
        state[j1] = 0
    end
end

function simulate!(state::Vector{Int}, N::Int, T::Int, p::Float64, alpha::Float64; periodic::Bool=false, rng=Random.default_rng())
    # TODO: precompute tables
    if periodic
        if iseven(N)
            rmax = div(N, 2)
        else
            rmax = div(N-1, 2)
        end
        r = 1:rmax
        weights = r.^(-alpha)    
        if iseven(N)
            weights[rmax] = weights[rmax]/2
        end
        probabilities = weights/sum(weights)
        prob, alias = make_alias_table(probabilities)
        x1s, x2s = sample_pair_periodic(N, prob, alias, T; rng=rng)
    else
        r = 1:N
        weights = (N .- r).* (r.^(-alpha))
        probabilities = weights/sum(weights)
        prob, alias = make_alias_table(probabilities)
        x1s, x2s = sample_pair(prob, alias, T; rng=rng)
    end
    fill!(state, 0) # reinitialize state
    for t = 2:T
        x = rand(rng)
        if x >= p
            s1 = x1s[t]
            s2 = x2s[t]
            measure_two!(s1, s2, state)
        else
            s = rand(rng, 1:N)
            measure_one!(s, state)
        end
    end
    return state
end

function simulate(N::Int, T::Int, p::Float64, alpha::Float64; periodic::Bool=false, rng=Random.default_rng())
    state = zeros(Int, N)
    simulate!(state, N, T, p, alpha; periodic=periodic, rng=rng)
    return state
end

function simulate_multiple(N::Int, T::Int, p::Float64, alpha::Float64, trials::Int; periodic::Bool=false, threaded::Bool=false)
    states = Matrix{Int}(undef, N, trials)
    if threaded
        Threads.@threads for t in 1:trials
            rng = Xoshiro(1234+t)
            state = simulate(N, T, p, alpha; periodic=periodic, rng=rng)
            states[:,t] = state
        end
    else
        for t in 1:trials
            rng = Xoshiro(1234+t)
            states[:,t] = simulate(N, T, p, alpha; periodic=periodic, rng=rng)
        end
    end
    return states
end

function average_entropy_profile(N::Int, T::Int, p::Float64, alpha::Float64, trials::Int; periodic::Bool=false)
    thread_profiles = [zeros(Int, N-1) for _ in 1:Threads.nthreads()]
    thread_states = [zeros(Int, N) for _ in 1:Threads.nthreads()]
    Threads.@threads for t in 1:trials
        rng = Xoshiro(1234+t) 
        state = thread_states[Threads.threadid()]
        state = simulate!(state, N, T, p, alpha; periodic=periodic, rng=rng)
        entropy_profile!(thread_profiles[Threads.threadid()], state)
    end
    entropy_profile = sum(thread_profiles) / trials
    return entropy_profile
end