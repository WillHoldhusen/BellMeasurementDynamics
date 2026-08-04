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

