export placepoints

function placepoints(tm::TopoMapStruct, edges::Array{Int64,2}, weights::Vector{Float64})

    tm.comps = IntDisjointSets(size(edges, 1) + 1)
    tm.compMap = Vector{Component}(undef, size(edges, 1) + 1)
    tm.verts = Vector{Vertex}(undef, size(edges, 1) + 1)


    for i = 1:length(tm.compMap)
        v = Vertex(Point(0, 0), i)
        tm.verts[i] = v
        c = Component([v.id], [v.p, v.p])
        tm.compMap[i] = c
    end
    order = sortedges(weights)

    for oi = 1:length(order)
        i = order[oi]
        p1 = edges[i, 1]
        p2 = edges[i, 2]

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
    c1_vert_len = length(c1.vertices)
    c2_vert_len = length(c2.vertices)
    merged = Component(Vector{Int64}(undef, c1_vert_len + c2_vert_len), [])
    for i = 1:c1_vert_len
        merged.vertices[i] = c1.vertices[i]
    end
    for i = 1:c2_vert_len
        merged.vertices[i+c1_vert_len] = c2.vertices[i]
    end

    if edge_len > 0.0
        # Transform vertices of the two components appropriately 
        t1::Transformation = alignhull(c1.hull, tm.verts[v1].p, true) # aligned w.r.t. top edge, so stays in the bottom
        transformcomponents(tm, c1, t1, 0.0)

        t2::Transformation = alignhull(c2.hull, tm.verts[v2].p, false) # aligned w.r.t. bottom edge, so offset should be added
        transformcomponents(tm, c2, t2, edge_len)

        # Compute the hull of merged component
        pts = Array{Point}(undef, length(c1.hull) + length(c2.hull))

        c1_length = length(c1.hull)
        c2_length = length(c2.hull)
        for i = 1:c1_length
            pts[i] = transform(c1.hull[i], t1, 0.0)
        end
        for i = 1:c2_length
            pts[i+c1_length] = transform(c2.hull[i], t2, edge_len)
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
end
