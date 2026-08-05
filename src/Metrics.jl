function bonds_across(state, cut)
    return count(x -> x > cut, state[1:cut])
end

function bonds_outside(state, start, finish)
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

function unpaired_sites(state)
    return count(iszero, state)
end

function count_rs(state)
    N = length(state)
    rs = zeros(Int, N)
    for (i, s) in enumerate(state)
        if s != 0
            r = abs(i-s)
            rs[r] += 1
        end
    end
    return rs
end

