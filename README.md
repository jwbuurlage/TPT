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
- (optional) [zeromq](zeromq.org) for communicating with visualization servers
- (optional) [boost::hana](http://www.boost.org/doc/libs/1_61_0/libs/hana/doc/html/index.html) is used to generate the Python bindings.
- (optional) [CUDA](http://www.nvidia.com/object/cuda_home_new.html) (>= 7.0)

Provided as submodules
- [Catch](https://github.com/philsquared/Catch), for unit tests
- [fmt](https://github.com/fmtlib/fmt) as an iostream replacement
- [bulk](https://github.com/jwbuurlage/Bulk) for distributed computing
- (optional) [slicevis](https://github.com/jwbuurlage/slicevis) as an visualization server
- (optional) [pybind11](https://github.com/pybind/pybind11) to generate Python bindings

The following build tools should be available:
- [CMake](https://cmake.org/) (>= 3.0)
- Modern C++ compiler (with support for at least C++14), e.g. GCC >= 6.0 or clang >= 3.8

The library is being tested on Fedora 24 and Arch Linux, but the code should be portable to other platforms.

# Building

## Examples

The core of the library is header only, so it does not have to be built itself. We start with initializing the submodules:

```
git submodule init
git submodule update --remote
```

To build the examples:

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

The Python bindings can be used through `tomo.py` in the `Python` folder, which also adds some rudimentary plotting functionality on top of the bindings.

# Building with optional features

To build with ZMQ (and MPI), run instead of `cmake .`:

```
cmake -DWITH_MPI=on -DWITH_ZMQ=on .
```

CUDA support has been disabled for now.

# Writing your own algorithms

## C++

First of all, have a look at the examples provided in the `examples` folder.

There are four core components that are used for reconstruction:

- `tomo::volume` represents the volume geometry, i.e. the part of space in which the image object resides.
- `tomo::dim` is the namespace for the 'discrete integration methods' (also called interpolators, projectors or kernels).
    - `closest` projects any point on the ray to the closest voxel
    - `linear` does D-dimensional linear interpolation around the ray point to the surrounding voxels
    - `joseph` does (D-1) dimensional linear interpolation by considering points on the ray that have integer coordinates in one fixed dimension.
- `tomo::geometry` is the namespace for the various acquisition geometries
    - `cone_beam`
    - `dual_axis_parallel`
    - `dynamic_cone_beam`
    - `fan`
    - `helical_cone_beam`
    - `laminography`
    - `list` is a list of lines without any implied stucture.
    - `parallel<2_D>`
    - `parallel<3_D>`
    - `tomosynthesis`
    - `trajectory` is the base class for cone-beam-like geometries where the source and the (position and tilt of the) detector follow a given path.
- `tomo::image` represents the image data, there is only one phantom
    - `shepp_logan<2_D>`
    - `shepp_logan<3_D>`
    - `modified_shepp_logan<2_D>`
    - `modified_shepp_logan<3_D>`

These can be used together completely independently. The reason is that we take the following approach to these concepts:

- A geometry acts as nothing more than a container of lines, so you can write:
```
for (auto line : geometry) {
    // use line
}
```
- A discrete integration method takes a line, and produces a number of 'matrix elements', that contain the voxel (as an index), and the attenuation coefficient (value of the matrix element):
```
for (auto element : projector(line)) {
    // element.index is the voxel
    // element.value is the coefficient
}
```
Using this approach, many interesting algorithms can be written in an efficient but flexible manner.

There are also some standard algorithms implemented, including `ART`, `SART`, and `SIRT`.

## Python

The Python bindings expose the different concepts (images, volumes, geometries and dims) as well as the standard implemented algorithms.
