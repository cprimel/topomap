using TopoMap
using Test, Distributions, MLDatasets, LazySets
using EMST


function test_project()
    data = rand(Uniform(0, 10), (100,4))

    tm = TopoMap.project(data)
    @test summary(tm) == "100×2 Matrix{Float64}"
    
    tm = TopoMap.project(data;dimensions=3,kernel_estimator="normal")
    @test summary(tm) == "100×3 Matrix{Float64}"

    #weights = Vector{Float64}([])
    #oldfromnew = Vector{Int64}([])


    #edges = Vector{Tuple{Int64,Int64}}(undef, size(e_out, 1))
    #print(size(e_out))
    #for i = 1:size(e_out, 1)
    #    indexA::Int64 = oldfromnew[e_out[i, :][1]]
    #    indexB::Int64 = oldfromnew[e_out[i, :][2]]

    #    if indexA < indexB
    #        edges[i] = (indexA, indexB)
    #    else
    #        edges[i] = (indexB, indexA)
    #    end

    #end


end


function test_find_angle(p1::Point, p2::Point)
    t = Transformation(0, 0, pi / 4, pi / 4)
    find_angle(p1, p2, t)
    @test t.cos == sqrt(1 - (pi / 4)^2)
    @test t.sin == -pi / 4
end

"""
Note: Sanity check not robust test!
"""
function test_compute_convex_hull()
    gen_pts = N -> [randn(2) for i in 1:N]
    v = gen_pts(30)

    points = Vector{Point}([])

    for i in eachindex(v)
        push!(points, Point(v[i][1], v[i][2]))
    end

    chull = Vector{Point}([])

    compute_convex_hull(points, chull)
    # Test to see if result includes points with min and max
    min_x = Inf
    min_y = Inf
    max_x = -Inf
    max_y = -Inf

    for p in points
        x_val = p.x
        y_val = p.y
        min_x = x_val < min_x ? x_val : min_x
        max_x = x_val > max_x ? x_val : max_x
        min_y = y_val < min_y ? y_val : min_y
        max_y = y_val > max_y ? y_val : max_y
    end

    contains_min_x = false
    contains_max_x = false
    contains_min_y = false
    contains_max_y = false

    for p in chull
        x_val = p.x
        y_val = p.y
        contains_min_x = (x_val == min_x) || contains_min_x
        contains_max_x = (x_val == max_x) || contains_max_x
        contains_min_y = (y_val == min_y) || contains_min_y
        contains_max_y = (y_val == max_y) || contains_max_y
    end
    @test contains_min_x && contains_max_x && contains_min_y && contains_max_y
end


function test_3d_project()
    data = rand(Uniform(0, 10), (4, 100))
    tm = TopoMap.project(data; dimensions=3)
    println(tm)
end

@testset "Utility tests" begin


    @testset "Test computeConvexHull" begin
        test_compute_convex_hull()
    end
    #@testset "Test IntDisjointSets" begin
    #    s = IntDisjointSets(10)
    #    x = find(s, 5)
    #    @test x == 5
    #    mergeset(s, 3, 5)
    #end

    @testset "Test project" begin
        test_project()
    end

    #@testset "Test alignhull" begin
    #   c1 = Component([1], [Point(1, 0), Point(1, 2), Point(0,0)])
    #  p1 = Point(0, 0)
    # println("Hull pre-align: $(c1.hull)")
    #trans = alignhull(c1.hull, p1, false)
    #println("Translation: $trans")
    #println("Hull post-align: $(c1.hull)")

    #for i = 1:length(c1.hull)
    #    c1.hull[i] = transform(c1.hull[i], trans, 0.0)
    #end
    #println("Hull post-transform: $(c1.hull)")

    #end
end