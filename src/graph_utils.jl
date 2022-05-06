export Vertex, Component, sortedges

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

function sortedges(weights::Array{Float64})
    order = collect(1:length(weights))
    perm = sortperm(weights)
    order = order[perm]
    return order
end