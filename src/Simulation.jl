include("Entropy.jl")
include("AliasTables.jl")


function measure_one!(s, state)
    j = state[s]
    if j != 0
        state[s] = 0
        state[j] = 0
    end
end


function measure_two!(s1, s2, state)
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


function simulate(N, T, p, alpha)
    r = 1:N
    weights = (N .- r).* (r.^(-alpha))
    
    probabilities = weights/sum(weights)

    prob, alias = make_alias_table(probabilities)
    x1s, x2s = sample_pair_open(prob, alias, T)
    state = zeros(Int, N)
    for t = 2:T
        x = rand()
        if x >= p
            s1 = x1s[t]
            s2 = x2s[t]
            measure_two!(s1, s2, state)
        else
            s = rand(1:N)
            measure_one!(s, state)
        end
    end
    return state
end

function simulate_multiple(N, T, p, alpha, trials)
    # entropies = zeros(Int, N, trials)
    # n1s = zeros(Int, threads)
    states = zeros(Int, trials, N)
    for t in 1:trials
        states[t,:] = simulate(N, T, p, alpha)
        # for i in range(N)
        #     entropies[i, t] = bonds_across(state, div(N, 2))
        # end
        # n1s[t] = unpaired_sites(state)
    end
    # entropy_avg = zeros(N)
    # for i in range(N)
    #     entropy_avg[i] = sum(entropies[i,:])/trials
    # end
    # n1_avg = sum(n1s)/trials
    # return entropy_avg, n1_avg
    return states
end

function simulate_multiple_threaded(N, T, p, alpha, trials)
    # entropies = zeros(Int, N, trials)
    # n1s = zeros(Int, threads)
    states = zeros(Int, trials, N)
    Threads.@threads for t in 1:trials
        states[t,:] = simulate(N, T, p, alpha)
        # for i in range(N)
        #     entropies[i, t] = bonds_across(state, div(N, 2))
        # end
        # n1s[t] = unpaired_sites(state)
    end
    # entropy_avg = zeros(N)
    # for i in range(N)
    #     entropy_avg[i] = sum(entropies[i,:])/trials
    # end
    # n1_avg = sum(n1s)/trials
    # return entropy_avg, n1_avg
    return states
end