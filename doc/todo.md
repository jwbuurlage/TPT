# CONTINUE:
* [x] Start implementing SART and SIRT
* [x] Add iterations to ART, how to add initial guess (which is now of type
    'std::optional') to Python (A: ignore initial guess for now)
* [x] Expand Python functionality (better plots, performance / convergence
    graphs)
* [ ] Make simple profiling environment
* [ ] CUDA / Distributed is next
* [ ] Optimize CPU algorithms (cache-use, limit recomputations, data-oriented, ...)
* [ ] 3D / cone / slices
* [ ] Future: work on visualizer? Parallel marching cubes? How to deal with massive amount of data on network
* [ ] We want the objects to be able to morph into matrices and vectors for algorithm operations
* [ ] Should cache 'R', 'C' and 'w_norms' somewhere so we can perform iterations independently


# List of ideas
* Vectors and intersection algorithms, kernels. Which file and what do we implement ourselves?
* Templated accelerator functions (or component like objects)
* Accelerator/platform agnostic
* Write also a separate 2D/3D visualizer (tomovis), agnostic renderer and can
    be off-site
* Use python for interacting, C++ for testing but should also be used standalone
    or as a library
* In order to accelerate efficiently, we need to layout object data efficiently
    (*data oriented development*) to prevent inefficient copying
* Do we want to distuingish between *volume* and *image*
* Simulate noise?
* Templated images for double/single precision? May conflict with CUDA / OpenCL.
* One 'geometry' that uses an iterator to generates line, one 'kernel iterator' that generates the triplets from lines

# Glossary of terms
* *Volume*
* *Image*
* *Geometry*
* *Line*, initial, delta
* *Line iterator*
* *Kernel iterator*

# Preparation and utilities
* [ ] Load or generate phantom images
* [ ] Geometry objects
    - These should be able to structurally generate lines
    - Should also be able to be an arbitrary collection of lines
* [ ] Think about how to incorporate distributed algorithms from the start
* [x] Need to interface with python so that we can plot and test, pybind11
* [ ] Run unit tests in C++, interact through python front-end

# Algorithms
* [ ] Need accelerated methods for forward and backward projection for arbitrary
        geometry and dimension
* [ ] Can we construct ARMs from these smaller primitives (scaling, forward,
        backward) or need custom
* [ ] Use OpenCL or still CUDA
