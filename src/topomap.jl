export TopoMapStruct

"""
Keeps track of connected components of graph and cumulative transformations.

...
# Arguments
- `leafSize::Integer`:
...
"""
mutable struct TopoMapStruct

    data::Array{Float64,2}
    dimensions::Int64
    leafSize::Int64
    verbose::Bool
    comps::IntDisjointSets
    compMap::Vector{Component}
    verts::Vector{Vertex}

    TopoMapStruct(data, dim, leafSize, verbose) = new(data, dim, leafSize, verbose)
end