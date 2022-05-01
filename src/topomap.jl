# using DataStructures
using EMST

export TopoMapStruct, project, placepoints, mergecomponents, transformcomponents, sortedges, transform

import Base.isequal

"""
    TopoMap(data, dim, leafSize, verbose)
 
Keeps track of connected components of graph.

...
# Arguments
- `leafSize::Integer`:
...
"""
mutable struct TopoMapStruct

    data::Array{Float64,2}
    dim::Int64
    leafSize::Int64
    verbose::Bool
    comps::IntDisjointSets
    compMap::Vector{Component}
    verts::Vector{Vertex}
    #steps::Vector{Vector{Vertex}}

    TopoMapStruct(data, dim, leafSize, verbose) = new(data, dim, leafSize, verbose)
end

function project(data::Array{Float64,2}, leafSize::Int64, verbose::Bool)
    res = Array{Float64,2}(undef, 2, size(data, 2))
    tm = TopoMapStruct(data, 2, leafSize, verbose)
    weights = Vector{Float64}([])
    oldfromnew = Vector{Int64}([])

    println("Computing EMST...")
    edges, weights, oldfromnew = compute_emst(tm.data; leafSize)

    println("Placing points...")

    pts = placepoints(tm, edges, weights)

    println(size(pts))
    println(size(res))
    for i in eachindex(pts)
        res[1, i] = pts[i].x
        res[2, i] = pts[i].y
    end
    return res#, tm.steps
end


function placepoints(tm::TopoMapStruct, edges::Array{Int64,2}, weights::Vector{Float64})

    tm.comps = IntDisjointSets(size(edges,1) + 1)
    tm.compMap = Vector{Component}(undef, size(edges,1) + 1)
    tm.verts = Vector{Vertex}(undef, size(edges,1) + 1)
    #tm.steps = Vector{Vector{Vertex}}()

    for i = 1:length(tm.compMap)
        v = Vertex(Point(0, 0), i)
        tm.verts[i] = v
        c = Component([v.id], [v.p, v.p])
        tm.compMap[i] = c
    end
    order = sortedges(weights)

    for oi = 1:length(order)
        i = order[oi]
        p1 = edges[i,1]
        p2 = edges[i,2]

        c1 = find(tm.comps, p1)
        c2 = find(tm.comps, p2)
        if c1 == c2
            println("Error! MST edge belongs to same component!")
            return
        end

        comp1::Component = tm.compMap[c1]
        comp2::Component = tm.compMap[c2]

        comp::Component = mergecomponents(tm, comp1, comp2, p1, p2, weights[i])
        merge(tm.comps, c1, c2)
        c = find(tm.comps, c1)
        tm.compMap[c] = comp
    end

    pts = Vector{Point}(undef, length(tm.verts))
    for i = 1:length(tm.verts)
        pts[i] = tm.verts[i].p
    end
    return pts
end

function mergecomponents(tm::TopoMapStruct, c1::Component, c2::Component, v1::Int64, v2::Int64, edge_len::Float64)
    # Compute merged component
    merged = Component([], [])
    append!(merged.vertices, c1.vertices)
    append!(merged.vertices, c2.vertices)

    if edge_len > 0.0
        # Transform vertices of the two components appropriately 
        t1::Transformation = alignhull(c1.hull, tm.verts[v1].p, true) # aligned w.r.t. top edge, so stays in the bottom
        transformcomponents(tm, c1, t1, 0.0)

        t2::Transformation = alignhull(c2.hull, tm.verts[v2].p, false) # aligned w.r.t. bottom edge, so offset should be added
        transformcomponents(tm, c2, t2, edge_len)

        # Compute the hull of merged component
        pts = Vector{Point}([])

        for i = 1:length(c1.hull)
            push!(pts, transform(c1.hull[i], t1, 0.0))
        end
        for i = 1:length(c2.hull)
            push!(pts, transform(c2.hull[i], t2, edge_len))
        end

        compute_convex_hull(pts, merged.hull)
    else
        if length(c1.hull) != 2 || length(c2.hull) != 2
            println("Error! Hull cannot have more than one point when edge length is 0.")
        end

        push!(merged.hull, c2.hull[1])
        push!(merged.hull, c2.hull[2])
    end
    return merged
end


function transformcomponents(tm::TopoMapStruct, c::Component, t::Transformation, y_offset::Float64)

    for i = 1:length(c.vertices)
        vin::Int64 = c.vertices[i]
        tm.verts[vin].p = transform(tm.verts[vin].p, t, y_offset)
    end
    #push!(tm.steps, deepcopy(tm.verts))
end
