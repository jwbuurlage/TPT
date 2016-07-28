# TODAY'S GOAL
[ ] Python interface to CUDA backend
[ ] Add CLI args to reconstruction cpp examples
[ ] Start optimizing CUDA SART

[x] Implement 'closest kernel' in CPU, see if we get the same result (ANS: not quite, but fixed a gpu bug)
[x] Clean up GPU code, think of simple abstractions

# CONTINUE:
* [ ] Set up unit tests, set up documentation
* [ ] Optimize CPU algorithms (cache-use, limit recomputations, data-oriented, ...)
* [ ] CUDA and Distributed is next
* [ ] Need *good* abstractions for GPU / Distributed objects (geometry, projectors, ...)
* [ ] 3D / cone / slices

# List of ideas
* We want the objects to be able to morph into matrices and vectors for algorithm operations
* Should cache 'R', 'C' and 'w_norms' somewhere so we can perform iterations independently for testing. This can be done by objectifying the reconstruction algorithms.
* Vectors and intersection algorithms, kernels. Which file and what do we implement ourselves?
* Templated accelerator functions (or component like objects)
* Accelerator/platform agnostic
* Write also a separate 2D/3D visualizer (tomovis), agnostic renderer and can be off-site, see TomViz by cornell. Parallel marching cube.
* Use python for interacting, C++ for testing but should also be used standalone or as a library
* In order to accelerate efficiently, we need to layout object data efficiently (*data oriented development*) to prevent inefficient copying
* Do we want to distuingish between *volume* and *image*
* Simulate noise?
* Templated images for double/single precision? May conflict with CUDA / OpenCL (can specialize for these systems).
* One 'geometry' that uses an iterator to generates line, one 'kernel iterator' that generates the triplets from lines.

# Preparation and utilities
* [ ] Run unit tests in C++, interact through python front-end

# Known bugs:
- SART fails (diverges) when using linear projector on over determined system.
