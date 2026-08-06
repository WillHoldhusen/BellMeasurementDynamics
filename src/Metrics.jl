function entropy_profile(state)
    N = length(state)
    entropies = zeros(Int, N-1)
    num_crossings = 0
    for i in 1:N-1
        j = state[i]
        if j > i
            # bond is crossing, add one
            num_crossings += 1
        elseif j != 0
            # second site of bond is to the left, no longer crossing 
            num_crossings -= 1
        end
        entropies[i] = num_crossings
    end
    return entropies
end


function entropy_profile!(old_profile, state)
    N = length(state)
    num_crossings = 0
    for i in 1:N-1
        j = state[i]
        if j > i
            # bond is crossing, add one
            num_crossings += 1
        elseif j != 0
            # second site of bond is to the left, no longer crossing 
            num_crossings -= 1
        end
        old_profile[i] += num_crossings
    end
    return old_profile
end


function bonds_across(state::Vector{Int}, cut::Int)
    return count(x -> x > cut, state[1:cut])
end

function bonds_outside(state::Vector{Int}, start::Int, finish::Int)
    bonds = 0
    for s in state[start:finish]
        if s > finish
            bonds += 1
        elseif s != 0 && s < start
            bonds += 1
        end
    end
    return bonds
end

function unpaired_sites(state::Vector{Int})
    return count(iszero, state)
end

function count_rs(state::Vector{Int})
    N = length(state)
    rs = zeros(Int, N-1)
    for i in 1:N
        j = state[i]
        if j > i
            rs[j-i] += 1
        end
    end
    return rs
end

function effective_exponents(Ns, SNs)
    return log.(SNs[2:end] ./ SNs[1:end-1]) ./ log.(Ns[2:end] ./ Ns[1:end-1])
end