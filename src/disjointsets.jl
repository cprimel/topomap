
export IntDisjointSets, merge, find, mergeset

mutable struct IntDisjointSets
    set::Array{Int64}
    function IntDisjointSets(size::Int64)
        s = new()
        s.set = fill(-1,1,size)
        return s
    end
end

function merge(s::IntDisjointSets, el1::Int64, el2::Int64)
    mergeset(s, find(s, el1), find(s, el2))
end

function find(s::IntDisjointSets, x::Int64)
    f::Int64 = s.set[x]
    if f < 0
        return x
    else
        xx::Int64 = find(s, f)
        s.set[x] = xx
        return xx
    end

end

function mergeset(s::IntDisjointSets, root1::Int64, root2::Int64)

    if root1 == root2
        return
    end
    r1::Int64 = s.set[root1]
    r2::Int64 = s.set[root2]
    if r2 < r1
        s.set[root1] = root2
    else
        if r1 == r2
            r1_stepback = r1 - 1
            s.set[root1] = r1_stepback
        end
        s.set[root2] = root1
    end
end