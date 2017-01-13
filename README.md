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
- [glm](http://glm.g-truc.net/0.9.8/index.html) header only mathematics library mimicking GLSL
- [boost::program_options](boost.org) for portable program options
- (optional) [boost::hana](http://www.boost.org/doc/libs/1_61_0/libs/hana/doc/html/index.html) is used to generate the Python bindings.
- (optional) [CUDA](http://www.nvidia.com/object/cuda_home_new.html) (>= 7.0)

Provided as submodules
- [Catch](https://github.com/philsquared/Catch), for unit tests
- [fmt](https://github.com/fmtlib/fmt) as an iostream replacement
- [bulk](https://github.com/jwbuurlage/Bulk) for distributed computing
- [slicevis](https://github.com/jwbuurlage/slicevis) as an optional visualization server
- (optional) [pybind11](https://github.com/pybind/pybind11) to generate Python bindings

The following build tools should be available:
- [CMake](https://cmake.org/) (>= 3.0)
- Modern C++ compiler (with support for at least C++14)

The library is being tested on Fedora 24 and Arch Linux, but the code should be portable to other platforms.

# Building

## Examples

The library is header only, so it does not have to be built itself. To build the examples:

```
cd build
cmake .
make
```

The resulting binaries will be in the `bin` folder.

# Building the Python bindings

To generate the Python bindings:

```
cd python/build
cmake .
make
```

The python bindings can be used through `tomo.py` in the `Python` folder, which also adds some rudimentary plotting functionality on top of the bindings.

# Writing your own algorithms

## C++

## Python
