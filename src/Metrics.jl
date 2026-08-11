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


function bonds_outside(state::Vector{Int}, lo::Int, hi::Int)
    bonds = 0
    for i in lo:hi
        s = state[i]
        if s != 0 && (s < lo || s > hi)
            bonds += 1
        end
    end
    return bonds
end

function is_outside(i, lo, hi)
    if i != 0 && !(lo <= i <= hi)
        return 1
    else
        return 0
    end
end

function central_entropy_profile!(entropy::Vector{Int}, state::Vector{Int})
    N = length(state)
    lo = div(N, 2)
    hi = lo + 1
    entropy[1] = is_outside(state[lo], lo, hi)
    entropy[1] += is_outside(state[hi], lo, hi)

    for k in 2:div(N,2)
        new_lo = lo - 1
        new_hi = hi + 1

        i = state[new_lo]
        j = state[new_hi]

        if i != 0
            was_inside = lo <= i <= hi
            is_inside = new_lo <= i <= new_hi
            entropy[k] = entropy[k-1] - was_inside + !is_inside
        else
            entropy[k] = entropy[k-1]
        end
        if j != 0
            was_inside = lo <= j <= hi
            is_inside = new_lo <= j <= new_hi
            entropy[k] += -was_inside + !is_inside
        end
        lo = new_lo
        hi = new_hi
    end
    return entropy
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