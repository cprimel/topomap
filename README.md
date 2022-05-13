# TopoMap.jl

## Instructions

Requires [Boost](https://www.boost.org/) and [mlpack](https://www.mlpack.org/) headers. 

TopoMap is accessed through a single function `project` which takes as its only required input a `Matrix{Float64}`. Optional parameters include `dimensions::Int64=2`, `leafSize::Int64=1`, `kernel_estimator::String="none"`, `bandwidth::Tuple{Float64, Float64}=(NaN,NaN)`, and `verbose::Bool=false`

For examples of how to use TopoMap, see `examples\examples.ipynb`.

## Todos

Performance:
* Rewrite API to split `project` from 3D projection routines.
* Minimizing allocations:
  - [ ] Move to StaticArrays for Points
  - [ ] Minimize use of dynamic allocation (unsized vectors, push, append, etc.)