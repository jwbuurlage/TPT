# TODAY'S GOAL
* [ ] Learn MPI, write Bulk backend

# CONTINUE:
* [ ] A 'CUDA stream' can have associated buffers (geometry, image, sino bridges), which stay alive while the stream stays alive.
* [ ] Set up unit tests, set up documentation
* [ ] Optimize CPU algorithms (cache-use, limit recomputations, data-oriented, ...)
* [ ] CUDA and Distributed is next
* [ ] 3D / cone / slices
* [x] Need *good* abstractions for GPU / Distributed objects (geometry, projectors, ...)

# List of ideas
* We want the objects to be able to 'morph' into matrices and vectors for algorithm operations
* Set up proper profiling code, look into `boost::accumulator` (see Benchmarking C++ @ CppCon 2015).
* Should cache 'R', 'C' and 'w_norms' somewhere so we can perform iterations independently for testing. This can be done for example by objectifying the reconstruction algorithms.
* Vectors and intersection algorithms, kernels. Which file and what do we implement ourselves?
* Templated accelerator functions (or component like objects)
* Accelerator/platform agnostic
* Write also a separate 2D/3D visualizer (tomovis), agnostic renderer and can be off-site, see TomViz by cornell. Parallel marching cube.
* In order to accelerate efficiently, we need to layout object data efficiently (*data oriented development*) to prevent inefficient copying
* Simulate noise?

# Preparation and utilities
* [ ] Run unit tests in C++, interact through python front-end

# Known bugs:
- SART fails (diverges) when using linear projector on over-determined system.
