# Short-term:

It is necessary for my experiments to have distributed 3d reconstruction of arbitrary geometry. The inbetween steps:
- [x] Fix border error
- [x] Add tomo::args to distributed reconstruction
- [x] Have a version working for e.g. parallel geometry in 3d for ShLo phantom
	- [x] Set up a benchmarking tool (maybe adaptation of the Zee one)
    - [x] Need to plot the resulting image (per iteration), could maybe use Cinder for this
    - [x] Need to heavily optimize the 'harmonizing' (probably MPI backend of Bulk)
    - [x] We require a more general notion of a distributed sinogram (= stack of projections) that can be harmonized (simply rename)
    - [x] Add 3D partitioned phantom support
- [x] Make Joseph for D > 2
- [x] Update partitioning algorithm for new Bulk system.
    - [ ] Optimize the steps
    - [ ] Perform on downsampled geometries
    - [ ] Add voxel weight cumulative sum support
- [ ] Make the partitioning persistent, define some TOML format
    - [x] Integrate https://github.com/jbeder/yaml-cpp
    - [x] Alternative: TOML https://github.com/skystrife/cpptoml
    - [ ] Define TOML format for geometries
    - [ ] A problem is: how to pass geometry parameters, it is potentially different for each geometry.
    - [!] Add 'factory' for projectors and geometry, maybe as part of tomo::args
- [x] Support physical volume dimensions
    - [x] Volume should be subpixels (i.e. origin and lengths are `T` instead of `int`)
    - [x] All relative geometries and math functions should accomodate this change
    - [x] Fix DIMs working on voxels now that lengths is physical
- [ ] Make also a scheme for storing performance and partitioning results
    - [ ] Partitioning: communication volume
    - [ ] Performance: actually running the code
        - [ ] Version number software
        - [ ] Identifier cluster
- [ ] Start optimization: Restrict the geometry to the local volume
- [ ] Fix center error, interpolation shenanigans. (One of) the problem(s) is that 'inside' is relative to the projector.
      [x] This means that even though the interopolation point can be outside the volume, one of the 'hit' voxels can be inside, depending on the projector.
      [x] Should see if entering/exiting happens properly
      [x] For some projectors, do -1 delta and +1 delta
      [x] What if line misses, but one of the points hits the subvolume
      [x] Fixed a number of issues, but still a slight boundary error somewhere
- [ ] Real data:
    - [ ] Read e.g. data-exchange format
    - [ ] Read TIFF stack
- [ ] Install GCC 7.0 and use C++17 to clean up code

# Future:

- [x] Move to common geometry, no longer needed as template argument (needed for 'real data')
- [x] Move to common projector, no longer needed as template argument
- [ ] We require realistic parameters for geometries, and test 3D sirt with them
- [ ] Documentation for conventions (volume, voxels, origin, geometries, 1d vs multi_indices, ...)
- [ ] General interpolation scheme for Simon
- [ ] May want to test Bulk cpp backend too
- [ ] Would be nice to add to the benchmarking a way to visualize it in slicevis
- [ ] Optimize CPU algorithms (cache-use, limit recomputations, data-oriented, SIMD, ...)

# GPU
- [ ] Start working on new GPU system
    - [ ] A 'CUDA stream' can have associated buffers (geometry-, image-, sino bridges), which stay alive while the stream stays alive.
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
