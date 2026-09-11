using Random
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

function simulate!(state::Vector{Int}, N::Int, T::Int, p::Float64, alpha::Float64; 
                   periodic::Bool=false, rng=Random.default_rng(),
                   save_states::Bool=false)
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
    fill!(state, 0) # reinitialize 
    if save_states
        states = zeros(Int, N, T)
        states[:,1] = copy(state)
    end
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
        if save_states
            states[:,t] = copy(state)
        end
    end
    if save_states
        return state, states
    else
        return state
    end
end

function simulate_time_series(N::Int, T::Int, p::Float64, alpha::Float64; 
                              periodic::Bool=false, rng=Random.default_rng())
    # TODO: precompute tables
    state = zeros(Int, N)
    N0 = zeros(Int, T)
    N0[1] = N
    S = zeros(Int, T)
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
        N0[t] = count(iszero, state)
        c = 0
        for i = 1:div(N,2)
            c += state[i]>div(N,2)
        end
        S[t] = c
    end
    return N0, S
end

function simulate_track(N::Int, T::Int, p::Float64, alpha::Float64; periodic::Bool=false, rng=Random.default_rng())
    state = zeros(Int, N)
    r_direct = zeros(Int, N)
    r_indirect = zeros(Int, N)
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
            r_direct[abs(s1-s2)] += 1
            if state[s1] != 0 && state[s2] != 0
                r_indirect[abs(state[s1]-state[s2])] += 1
            end
            measure_two!(s1, s2, state)
        else
            s = rand(rng, 1:N)
            measure_one!(s, state)
        end
    end
    return state, r_direct, r_indirect
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

function simulate_multiple_time_series(N::Int, T::Int, p::Float64, alpha::Float64, trials::Int; periodic::Bool=false)
    N0s_thread = [zeros(Int, T) for _ in 1:Threads.nthreads()]
    Ss_thread = [zeros(Int, T) for _ in 1:Threads.nthreads()]
    Threads.@threads for t in 1:trials
        tid = Threads.threadid()
        rng = Xoshiro(1234+t)
        N0, S = simulate_time_series(N, T, p, alpha; periodic=periodic, rng=rng)
        N0s_thread[tid] .+= N0
        Ss_thread[tid] .+= S
    end
    N0s = sum(N0s_thread) ./ trials
    Ss = sum(Ss_thread) ./ trials
    return N0s, Ss
end

function average_entropy_profile(N::Int, T::Int, p::Float64, alpha::Float64, trials::Int; periodic::Bool=false, centered::Bool=false)
    if centered
        thread_profiles = [zeros(Int, div(N,2)) for _ in 1:Threads.nthreads()]
    else
        thread_profiles = [zeros(Int, N-1) for _ in 1:Threads.nthreads()]
    end

    thread_states = [zeros(Int, N) for _ in 1:Threads.nthreads()]
    Threads.@threads for t in 1:trials
        rng = Xoshiro(1234+t) 
        state = thread_states[Threads.threadid()]
        state = simulate!(state, N, T, p, alpha; periodic=periodic, rng=rng)
        if centered
            entropy = zeros(Int, div(N,2))
            central_entropy_profile!(entropy, state)
            thread_profiles[Threads.threadid()] .+= entropy
        else
            entropy = add_entropy_profile!(thread_profiles[Threads.threadid()], state)
        end
    end
    entropy_profile = sum(thread_profiles)/trials
    
    return entropy_profile
end

function get_entropies_and_stats(N::Int, T::Int, p::Float64, alpha::Float64, trials::Int; periodic::Bool=false)
    thread_profiles = [zeros(Int, N-1) for _ in 1:Threads.nthreads()]
    thread_entropy  = [zeros(Int, N-1) for _ in 1:Threads.nthreads()]
    thread_S = [0 for _ in 1:Threads.nthreads()]
    thread_S2 = [0 for _ in 1:Threads.nthreads()]

    thread_N0 = [0 for _ in 1:Threads.nthreads()]
    thread_N02 = [0 for _ in 1:Threads.nthreads()]

    thread_states = [zeros(Int, N) for _ in 1:Threads.nthreads()]

    thread_rs = [zeros(Int, N-1) for _ in 1:Threads.nthreads()]

    Threads.@threads for t in 1:trials
        tid = Threads.threadid()
        #rng = Xoshiro(1234+t) 
        rng = Xoshiro(12345678 + t)
        state = thread_states[tid]
        state = simulate!(state, N, T, p, alpha; periodic=periodic, rng=rng)
        entropy = new_entropy_profile!(thread_entropy[tid], state)
        thread_profiles[tid] .+= entropy
        thread_S[tid] += entropy[div(N,2)]
        thread_S2[tid] += entropy[div(N,2)]^2

        N0 = count(iszero, state)
        thread_N0[tid] += N0
        thread_N02[tid] += N0^2

        rs = count_rs(state)
        thread_rs[tid] += rs
    end
    entropy_profile = sum(thread_profiles)/trials
    S_sum  = sum(thread_S)
    S2_sum = sum(thread_S2)

    S_mean = S_sum / trials
    S_var  = (S2_sum - trials * S_mean^2) / (trials - 1)
    S_sem  = sqrt(S_var / trials)

    N0_sum = sum(thread_N0)
    N02_sum = sum(thread_N02)

    N0_mean = N0_sum / trials
    N0_var = (N02_sum - trials * N0_mean^2) / (trials - 1)
    N0_sem = sqrt(N0_var / trials)

    rs_mean = sum(thread_rs)/trials
    return entropy_profile, S_mean, S_sem, N0_mean, N0_sem, rs_mean
end


function average_r_dist(N::Int, T::Int, p::Float64, alpha::Float64, trials::Int; periodic::Bool=false)
    thread_rs = [zeros(Int, N-1) for _ in 1:Threads.nthreads()]
    thread_states = [zeros(Int, N) for _ in 1:Threads.nthreads()]
    Threads.@threads for t in 1:trials
        tid = Threads.threadid()
        rng = Xoshiro(1234+t) 
        state = thread_states[tid]
        state = simulate!(state, N, T, p, alpha; periodic=periodic, rng=rng)
        rs = count_rs(state)
        thread_rs[tid] += rs
    end
    return sum(thread_rs)/trials
end