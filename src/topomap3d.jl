using Polyhedra

export placepoints3d

function placepoints3d(tm::TopoMapStruct, edges::Vector{Tuple{Int64,Int64}}, weights::Vector{Float64})

    tm.comps = IntDisjointSets(length(edges) + 1)
    tm.compMap = Vector{Component}(undef, length(edges) + 1)
    tm.verts = Vector{Vertex}(undef, length(edges) + 1)
    tm.steps = Vector{Vector{Vertex}}()

    for i = 1:length(tm.compMap)
        v = Vertex(Point(0, 0, 0), i)
        tm.verts[i] = v
        c = Component([v.id], [v.p, v.p])
        tm.compMap[i] = c
    end
    order = sortedges(weights)

    for oi = 1:length(order)
        i = order[oi]
        p1 = edges[i][1]
        p2 = edges[i][2]

        c1 = find(tm.comps, p1)
        c2 = find(tm.comps, p2)
        if c1 == c2
            println("Error! MST edge belongs to same component!")
            return
        end

        comp1::Component = tm.compMap[c1]
        comp2::Component = tm.compMap[c2]

        comp::Component = mergecomponents3d(tm, comp1, comp2, p1, p2, weights[i])
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

function mergecomponents3d(tm::TopoMapStruct, c1::Component, c2::Component, v1::Int64, v2::Int64, edge_len::Float64)
    # Compute merged component
    merged = Component([], [])
    append!(merged.vertices, c1.vertices)
    append!(merged.vertices, c2.vertices)

    if edge_len > 0.0
        # Transform vertices of the two components appropriately 
        t1::Transformation = alignhull3d(c1.hull, tm.verts[v1].p, true) # aligned w.r.t. top edge, so stays in the bottom
        transformcomponents(tm, c1, t1, 0.0)

        t2::Transformation = alignhull3d(c2.hull, tm.verts[v2].p, false) # aligned w.r.t. bottom edge, so offset should be added
        transformcomponents(tm, c2, t2, edge_len)

        # Compute the hull of merged component
        pts = Vector{Point}([])

        for i = 1:length(c1.hull)
            push!(pts, transform(c1.hull[i], t1, 0.0))
        end
        for i = 1:length(c2.hull)
            push!(pts, transform(c2.hull[i], t2, edge_len))
        end

        compute_convex_hull_3d(pts, merged.hull)
    else
        if length(c1.hull) != 2 || length(c2.hull) != 2
            println("Error! Hull cannot have more than one point when edge length is 0.")
        end

        push!(merged.hull, c2.hull[1])
        push!(merged.hull, c2.hull[2])
    end
    return merged
end

"""
3D Convex hull functions
"""
function alignhull3d(hull::Vector{Point}, p::Point, topEdge::Bool)
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
    v1 = Point(NaN, NaN, NaN)
    v2 = Point(NaN, NaN, NaN)
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
    trans = Transformation(NaN, NaN, NaN, NaN, NaN, NaN)
    trans.tx = -hull[v].x
    trans.ty = -hull[v].y
    trans.tz = -hull[v].z
    if length(hull) > 2
        find_angle(v1, v2, trans)
        trans.β = π / 4
    else
        trans.sin = 0
        trans.cos = 1
    end

    return trans
end


function compute_convex_hull_3d(pts::Vector{Point}, chull::Vector{Point})
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



    mpts = Matrix{Float64}(undef, 3, n)
    for i in eachindex(pts)
        mpts[1, i] = pts[i].x
        mpts[2, i] = pts[i].y
        mpts[3, i] = pts[i].z
    end
    #println("input: $mpts")
    #println("summary: $(summary(mpts))")
    v_x = vrep(transpose(mpts))
    removevredundancy(v_x, GLPK.Optimizer)
    println("$(summary(v_x))")
    hull = points(v_x)
    println("output: $hull.")
    println("summary: $(summary(hull))")
    #push!(hull, first(hull))
    #reverse!(hull)
    for x in hull
        push!(chull, Point(x[1], x[2], x[3]))
    end
end


"""
function test_compute_convex_hull_3d()
    x = rand(Float64, (3,100))

    points = Vector{Point}([])

    for i in 1:size(x,2)
        push!(points, Point(x[1,i], x[2,i], x[3,i]))
    end

    chull = Vector{Point}([])

    compute_convex_hull_3d(points, chull)

    # Test to see if result includes points with min and max
    min_x = Inf
    min_y = Inf
    min_z = Inf
    max_x = -Inf
    max_y = -Inf
    max_z = -Inf

    for p in points
        x_val = p.x
        y_val = p.y
        z_val = p.z
        min_x = x_val < min_x ? x_val : min_x
        max_x = x_val > max_x ? x_val : max_x
        min_y = y_val < min_y ? y_val : min_y
        max_y = y_val > max_y ? y_val : max_y
        min_z = z_val < min_z ? z_val : min_z
        max_z = z_val > max_z ? z_val : max_z
    end

    contains_min_x = false
    contains_max_x = false
    contains_min_y = false
    contains_max_y = false
    contains_min_z = false
    contains_max_z = false


    for p in chull
        x_val = p.x
        y_val = p.y
        z_val = p.z
        contains_min_x = (x_val == min_x) || contains_min_x
        contains_max_x = (x_val == max_x) || contains_max_x
        contains_min_y = (y_val == min_y) || contains_min_y
        contains_max_y = (y_val == max_y) || contains_max_y
        contains_min_z = (z_val == min_z) || contains_min_z
        contains_max_z = (z_val == max_z) || contains_max_z
    end
    @test contains_min_x && contains_max_x && contains_min_y && contains_max_y && contains_min_z && contains_max_z
end
"""