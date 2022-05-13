using KernelDensity
using mlpack: emst
#using EMST

export project

import Base.isequal


function log(message::String,verbose::Bool)
    if verbose
        println("[ Info ] $message")
    end
end

function project(data::Matrix{Float64}; dimensions::Int64=2, leafSize::Int64=1, kernel_estimator::String="none", bandwidth::Tuple{Float64,Float64}=(NaN, NaN), verbose::Bool=false)
    if kernel_estimator != "none" && dimensions != 3
        println("Incorrect parameters: Kernel density estimation is only available for 2-dimensional projection.")
    end


    tm = TopoMapStruct(data, dimensions, leafSize, verbose)
    res = Matrix{Float64}(undef, size(tm.data, 1), dimensions)

    log("Computing EMST...", verbose)
    edges = emst(data; leaf_size=leafSize, verbose=verbose)
    weights = edges[:,3]
    edges = trunc.(Int, edges[:,1:2]) .+1
    log("Placing points...", verbose)
    pts = placepoints(tm, edges, weights)
    for i in eachindex(pts)
        res[i, 1] = pts[i].x
        res[i, 2] = pts[i].y
    end
    if dimensions == 3
        if kernel_estimator == "normal"
            log("Estimating density...", verbose)
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
            res[1, 3] = (max_weight - weights[1]) / max_weight
            for i = 2:(size(res, 1)-2)
                res[i, 3] = (max_weight - (weights[i] + weights[i+1] / 2.0)) / max_weight
            end
            res[size(res, 1), 3] = (max_weight - weights[size(res, 1)-1]) / max_weight
        end

    end


    return res
end