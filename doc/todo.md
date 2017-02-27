# Short-term:

It is necessary for my experiments to have distributed 3d reconstruction given arbitrary geometry. The inbetween steps:

- [ ] Make a nice API that takes a data set, and reconstructs it in parallel
- [ ] Get a collection of realistic acquisition parameters
    - [ ] 'Update geometry format' for these parameters (i.e. just construct the files)
    - [ ] Update parser for these systems
    - [x] Update internal geometries for these systems
    - [ ] (Optional) update the Cinder based geometry / partitioning visualizer
- [ ] Update partitioning system
    - [ ] Redo trivial partitioning
        - [x] Update Bulk block partitionings to be able to choose the axes
        - [x] Use new Bulk partitioning
        - [ ] Function that computes overlap
        - [x] Choose best axis
    - [!] Update recursive bisectioning to new coordinate system
    - [ ] Make a scheme for partitioning results, add communication volume there
    - [ ] Perform the partitioning method on downsampled geometries
    - [ ] Run some partitioning examples
    - [ ] Allow to load the partitioning results and run with that, maybe as part of metadata
- [ ] Extend parallel beam definition, extend support in Tomos
    - [ ] Make a base class for trajectory-parallel beam for single- and dual-axis
- [ ] => RESULT: TABLE WITH COMMUNICATION VOLUMES VS SLABS (AND BLOCKS?)
- [ ] Run distributed performance tests
    - [ ] Make a scheme for results
        - [ ] Version number software
        - [ ] Identifier cluster
    - [ ] LaTeX tables and results, design experiments
- [ ] Performance optimization:
    - [ ] Restrict the geometry to the local volume
    - [ ] Extensive benchmarking, computation vs communication, timeline, etc.
- [ ] => RESULT: RUNTIMES AS FUNCTION OF COMMUNICATION VOLUME
- [ ] Optimize partitioning algorithm
    - [ ] Optimize the steps
    - [ ] Add voxel weight cumulative sum support
    - [ ] Add support for non-2^k-way partitionings
- [ ] => RESULT: ACTUAL TIMES ON CLUSTER (VOXEL?)
- [ ] Literature study
- [ ] Reading tomography papers
- [ ] Draft paper

# Future:

- [ ] Need some way for better debugging, visualization and inspection (e.g. detector positions to verify geometry etc.)
- [ ] Install GCC 7.0 and use C++17 to clean up code
- [ ] Bug: 3d plotter axis seem off (for lego?)
    - [ ] Documentation for internal representation conventions (volume, voxels, origin, geometries, 1d vs multi_indices, ...)
    - [ ] Verify nothing is mirrored or rotated
- [ ] We require realistic parameters for geometries, and test 3D sirt with them
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

# Done

- [x] Fix border error
- [x] Add tomo::args to distributed reconstruction
- [x] Have a version working for e.g. parallel geometry in 3d for ShLo phantom
	- [x] Set up a benchmarking tool (maybe adaptation of the Zee one)
    - [x] Need to plot the resulting image (per iteration), could maybe use Cinder for this
    - [x] Need to heavily optimize the 'harmonizing' (probably MPI backend of Bulk)
    - [x] We require a more general notion of a distributed sinogram (= stack of projections) that can be harmonized (simply rename)
    - [x] Add 3D partitioned phantom support
- [x] Make Joseph for D > 2
- [x] Fix center error, interpolation shenanigans. (One of) the problem(s) is that 'inside' is relative to the projector.
      [x] This means that even though the interopolation point can be outside the volume, one of the 'hit' voxels can be inside, depending on the projector.
      [x] Should see if entering/exiting happens properly
      [x] For some projectors, do -1 delta and +1 delta
      [x] What if line misses, but one of the points hits the subvolume
      [x] Fixed a number of issues, but still a slight boundary error somewhere
- [x] Real data:
    - [x] Read TIFF stack
- [x] Support physical volume dimensions
    - [x] Volume should be subpixels (i.e. origin and lengths are `T` instead of `int`)
    - [x] All relative geometries and math functions should accomodate this change
    - [x] Fix DIMs working on voxels now that lengths is physical
- [x] Run for (downsampled) medipix data, see if reconstruction works
    - [x] Add support for 'downsampling a reconstruction problem'
        - [x] Add changing voxels of a volume
        - [x] Add support for downsampling geometries, changing the detector pixels
        - [x] Downscale the projections to the new detector pixels
- [x] Move to common geometry, no longer needed as template argument (needed for 'real data')
- [x] Move to common projector, no longer needed as template argument
- [x] Update partitioning algorithm for new Bulk system.
- [x] Add geometry information to data, define some TOML format
    - [x] Alternative: TOML https://github.com/skystrife/cpptoml
    - [x] Define TOML format for geometries
        - [x] Add 'factory' for projectors and geometry, maybe as part of tomo::args
        - [x] Define cone beam definition, extend support in Tomos
            - [x] Redo trajectory geometry parameters
            - [x] Update all trajectories to new parent conventions
        - [x] Implement OpenCV TIFF reader
        - [x] Add support for loading projection stack from TIFF files to act as projections
            - [x] Make a 'reconstruction problem' struct

