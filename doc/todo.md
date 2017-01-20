# Short-term:

It is necessary for my experiments to have distributed 3d reconstruction of arbitrary geometry. The inbetween steps:
- [x] Fix border error
- [ ] Fix center error
- [ ] Have a version working for e.g. parallel geometry in 3d for ShLo phantom
	- [x] Set up a benchmarking tool (maybe adaptation of the Zee one)
    - [x] Need to plot the resulting image (per iteration), could maybe use Cinder for this
    - [x] Need to heavily optimize the 'harmonizing' (probably MPI backend of Bulk)
    - [ ] We require a more general notion of a distributed sinogram (!= stack of projections) that can be harmonized
    - [ ] Restrict the geometry to the local volume
- [ ] Read TIFF stack
- [ ] Read data-exchange format
- [ ] Install GCC 7.0 and use C++17 to clean up code
- [ ] Make Joseph for D > 2

# Future:
- [ ] Test bulk cpp backend too
- [ ] Optimize CPU algorithms (cache-use, limit recomputations, data-oriented, SIMD, ...)
- [ ] We may want a 'geometry::persist' that stores each geometry as a list geometry. This is more efficient then recomputing. This also makes distribution easier, because you can save the 'local origin' in this cache.
- [ ] Start working on GPU system
    - [ ] A 'CUDA stream' can have associated buffers (geometry, image, sino bridges), which stay alive while the stream stays alive.
    - [ ] CUDA storage backend
- [ ] Support for boxing and so on
- [ ] Non-local (boundary) operations

# List of ideas
* Set up proper profiling code, look into `boost::accumulator` (see Benchmarking C++ @ CppCon 2015).
* Should cache 'R', 'C' and 'w_norms' somewhere so we can perform iterations independently for testing. This can be done for example by objectifying the reconstruction algorithms.
* Vectors and intersection algorithms, kernels. Which file and what do we implement ourselves?
* Templated accelerator functions (or component like objects)
* Accelerator/platform agnostic
* In order to accelerate efficiently, we need to layout object data efficiently (*data oriented development*) to prevent inefficient copying
* Simulate noise

# Documentation
* For each entity (class, function, ...):
  - Brief description
  - (optional) detailed description, here high-level overview can be given (or intent).
  - (if relevant) Description of template parameters.
  - (if relevant) Description of parameters.
  - (if relevant) Description of return value.

# Known bugs:
- SART fails (diverges) when using linear projector on over-determined system.
