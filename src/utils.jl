using StaticArrays
using LazySets

export Point, Component, Vertex, Transformation, len, l2_squared, find_angle, compute_convex_hull, sortedges, alignhull

"""
Geometric utilities
"""

"""
Point
"""
mutable struct Point
    x::Float64
    y::Float64
    z::Float64
    Point(x, y) = new(x, y, 0.0)
    Point(x, y, z) = new(x, y, z)
end

function len(p::Point)::Float64
    sqrt(p.x^2 + p.y^2 + p.z^2)
end

function l2_squared(p1::Point, p2::Point)::Float64
    return (p1.x - p2.x)^2 + (p1.y - p2.y)^2 + (p1.z - p2.z)^2
end


"""
Vertex
"""
mutable struct Vertex
    p::Point
    id::Int64
end



"""
Component
"""
mutable struct Component
    vertices::Vector{Int64}
    hull::Vector{Point}
end

"""
Geometric transformations
"""

mutable struct Transformation
    tx::Float64
    ty::Float64
    cos::Float64
    sin::Float64
end


function transform(p::Point, t::Transformation, yOffset::Float64)
    x::Float64 = p.x + t.tx
    y::Float64 = p.y + t.ty

    xx::Float64 = x * cos(t.cos) - y * sin(t.sin)
    yy::Float64 = x * sin(t.sin) + y * cos(t.cos)
    
    yy += yOffset

    return Point(xx, yy)
end


function find_angle(p1::Point, p2::Point, t::Transformation)
    vec = Point(p2.x - p1.x, p2.y - p1.y)
    n::Float64 = len(vec)
    vec.x /= n
    vec.y /= n
    t.cos = vec.x
    t.sin = sqrt(1 - t.cos^2)
    if (vec.y) >= 0
        t.sin = -t.sin
    end
end

"""
2D Convex hull functions
"""
function alignhull(hull::Vector{Point}, p::Point, topEdge::Bool)
    # Find point in hull such that max_dist(point, p)
    v::Int64 = -1
    d2::Float64 = NaN
    for i = 1:length(hull)-1
        d::Float64 = l2_squared(hull[i], p)
        if v == -1
            d2 = d
            v = i
        else
            if d2 > d
                d2 = d
                v = i
            end
        end
    end
    # Note: Hull is ordered clockwise.
    v1 = Point(NaN, NaN)
    v2 = Point(NaN, NaN)
    if topEdge
        # make v, v+1 the top edge of the hull, s.t. v = (0,0)
        v1 = hull[v]
        v2 = hull[v+1]
    else
        # make v, v-1 the bottom edge of the hull, s.t. v = (0,0)
        if v == 1
            v = length(hull)
        end
        v1 = hull[v]
        v2 = hull[v-1]
    end
    trans = Transformation(NaN, NaN, NaN, NaN)
    trans.tx = -hull[v].x
    trans.ty = -hull[v].y
    if length(hull) > 2
        find_angle(v1, v2, trans)
    else
        trans.sin = 0
        trans.cos = 1
    end
    return trans
end


function compute_convex_hull(pts::Vector{Point}, chull::Vector{Point})
    empty!(chull)
    n = length(pts)

    if n == 1
        push!(chull, pts[1])
        push!(chull, pts[1])
        return
    elseif n == 2
        push!(chull, pts[1])
        push!(chull, pts[2])
        push!(chull, pts[1])
        return
    end

    mpts = Vector{Vector{Float64}}(undef, n)
    for i in eachindex(pts)
        mpts[i] = Vector{Float64}([pts[i].x, pts[i].y])
    end

    hull = convex_hull(mpts)

    # LazySets
    push!(hull, first(hull))
    reverse!(hull)
    for i in eachindex(hull)
        push!(chull, Point(hull[i][1], hull[i][2]))
    end
end



"""
Graph utilities
"""
function sortedges(weights::Array{Float64})
    order = collect(1:length(weights))
    perm = sortperm(weights)
    order = order[perm]
    return order
end