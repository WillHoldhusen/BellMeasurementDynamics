using Random

function make_alias_table(probabilities)
    N = length(probabilities)
    scaled = probabilities * N
    prob = zeros(N)
    alias = zeros(Int, N)
    small = []
    large = []
    for (i, p) in enumerate(scaled)
        if p < 1.0
            append!(small, i)
        else
            append!(large, i)
        end
    end
    while (!isempty(small) && !isempty(large))
        s = pop!(small)
        l = pop!(large)
        prob[s] = scaled[s]
        alias[s] = l

        scaled[l] = scaled[l] - (1.0 - scaled[s])

        if scaled[l] < 1.0
            append!(small, l)
        else
            append!(large, l)
        end
    end
    for i in vcat(large, small)
        prob[i] = 1.0
    end
    return prob, alias
end

function sample(prob, alias, size)
    N = length(prob)
    i = rand(1:N, size) # random integers between 1 and N 
    u = rand(size) # random floats in [0,1)
    return ifelse.(u .< prob[i], i, alias[i])
end

function sample_pair(prob, alias, size)
    N = length(prob)
    r = sample(prob, alias, size)
    x = [rand(1:N-d) for d in r]
    return x, x+r # pair order doesn't matter
end
