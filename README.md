# Tomos

Tomos is a library with many tools for tomographic reconstruction. Its features, and advantages over other 'tomography toolboxes', include:

- Modular system, where projectors, geometries, volumes and algorithmic skeletons are all independent and can be mixed and matched.
- Tomos supports a general model for distributed computing, enabling the resulting algorithms to run on clusters
- It is not bound to CUDA, MATLAB, or Python, although Python bindings and GPU computing are supported
- Written in modern C++, with OO / generic programming style like STL, enabling familiar abstractions without performance penalties.
- It has efficient CPU implementations for the algorithms.
- Support for 2D, 3D or higher dimensional reconstructions.

# Dependencies

The following libraries are required:

External:
- glm
- boost::program_options

Provided as submodules
- [Catch](https://github.com/philsquared/Catch), for unit tests
- [fmt](https://github.com/fmtlib/fmt) as an iostream replacement
- [bulk](https://github.com/jwbuurlage/Bulk) for distributed computing

For the Python bindings:
- pybind11
- boost::hana

To build:
- CMake (>= 3.0)
- Modern compiler (with support for at least C++14)

For GPU support:
- CUDA (>= 7.0)
