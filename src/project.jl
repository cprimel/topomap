using KernelDensity
using EMST

export project

import Base.isequal


function project(data::Matrix{Float64}; dimensions::Int64=2, leafSize::Int64=1, kernel_estimator::String="none", bandwidth::Tuple{Float64,Float64}=(NaN, NaN), verbose::Bool=false)
    if kernel_estimator != "none" && dimensions != 3
        println("Kernel density estimation is only available for 2-dimensional projection.")
    end


    tm = TopoMapStruct(data, dimensions, leafSize, verbose)

    res = Matrix{Float64}(undef, size(tm.data, 2), dimensions)
    weights = Vector{Float64}([])
    oldfromnew = Vector{Int64}([])

    println("Computing EMST...")
    edges, weights, oldfromnew = compute_emst(tm.data; tm.leafSize)

    println("Placing points...")
    if dimensions == 2
        pts = placepoints(tm, edges, weights)
        for i in eachindex(pts)
            res[i, 1] = pts[i].x
            res[i, 2] = pts[i].y
        end
    elseif dimensions == 3
        pts = placepoints(tm, edges, weights)
        for i in eachindex(pts)
            res[i, 1] = pts[i].x
            res[i, 2] = pts[i].y
        end
        if kernel_estimator == "normal"
            println("Estimating density...")
            if ~isnan(bandwidth[1]) && ~isnan(bandwidth[2])
                bi_dist = kde(res[:, 1:2]; bandwidth=bandwidth)
            else
                bi_dist = kde(res[:, 1:2])
            end
            for i = 1:size(res, 1)
                bd_x = searchsortednearest(collect(bi_dist.x), res[i, 1])
                bd_y = searchsortednearest(collect(bi_dist.y), res[i, 2])
                res[i, 3] = bi_dist.density[bd_x, bd_y]
            end
        else
            # Set Z to average weight on vertices edges
            max_weight = maximum(weights)
            sort!(weights)
            res[1, 3] = (max_weight - weights[1]) / max_weight
            for i = 2:(size(res, 1)-2)
                res[i, 3] = (max_weight - (weights[i] + weights[i+1] / 2.0)) / max_weight
            end
            res[size(res, 1), 3] = (max_weight - weights[size(res, 1)-1]) / max_weight
        end

    end


    return res
end