# TopoMap.jl

## Todos

Functionality
- [ ] 2d projection works and shows proper clustering
    - [x] 3 x spheres
    - [x] Iris 
    - [ ] more artificial and toy datasets to test on
- [ ] 3d projection
  - [x] convex hull computation set up (Polyhedra.jl)
  - [ ] integrate 3d projection into current code (first experiment showed no performance hit for branching)
  - [ ] test placement strategies   

Performance:
- [ ] Move to plain arrays or even StaticArrays.jl where applicable; a lot of time is spent allocating memory - expect that profiling would reveal functions with heavy data copying / pushing to vectors would be the primary culprit
