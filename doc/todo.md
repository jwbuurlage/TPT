#CONTINUE:
* [ ] A 'CUDA stream' can have associated buffers (geometry, image, sino bridges), which stay alive while the stream stays alive.
* [ ] Other geometries
    - [ ] 3D
    - [ ] fan
    - [x] slices (for parallel)
* [ ] Optimize CPU algorithms (cache-use, limit recomputations, data-oriented, SIMD, ...)
* [ ] We may want a 'geometry::persist' that stores each geometry as a list geometry. This is more efficient then recomputing. This also makes distribution easier, because you can save the 'local origin' in this cache.

#List of ideas
* We want the objects to be able to 'morph' into matrices and vectors for algorithm operations
* Set up proper profiling code, look into `boost::accumulator` (see Benchmarking C++ @ CppCon 2015).
* Should cache 'R', 'C' and 'w_norms' somewhere so we can perform iterations independently for testing. This can be done for example by objectifying the reconstruction algorithms.
* Vectors and intersection algorithms, kernels. Which file and what do we implement ourselves?
* Templated accelerator functions (or component like objects)
* Accelerator/platform agnostic
* We need also a separate 2D/3D visualizer (I once started working on tomovis), agnostic renderer and can be off-site, see TomViz by cornell. Parallel marching cube. This could be part of the 'tomo pipeline' idea.
* In order to accelerate efficiently, we need to layout object data efficiently (*data oriented development*) to prevent inefficient copying
* Simulate noise?

#Documentation
* For each entity (class, function, ...):
  - Brief description
  - (optional) detailed description, here high-level overview can be given (or intent).
  - (if relevant) Description of template parameters.
  - (if relevant) Description of parameters.
  - (if relevant) Description of return value.

#Known bugs:
- SART fails (diverges) when using linear projector on over-determined system.
