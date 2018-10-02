# Tomos

![](doc/images/logo.png)

Tomos is a software library and tools for tomographic reconstruction.

It uses a modular system, where kernels, geometries, algorithmic skeletons are independent and can be mixed and matched.
We support a general model for distributed computing, enabling the resulting algorithms to run on clusters
There is support for 2D, 3D or higher dimensional reconstructions, and support for arbitrary floating point types.

# Tools

- `tomos_partition` a partitioning tool based on the _geometric recursive coordinate bisectioning_ algorithm.

   It is also possible to perform partitioning on geometries from the ASTRA toolbox.
   This can be done using the Python bindings, see `python/examples/astra_partition.py`.

# Library

## Code example

Reconstructing a Shepp-Logan phantom using SIRT:

```cpp
#include "tomos/tomos.hpp"

int main() {
    using T = float;
    constexpr tomo::dimension D = 2_D;

    int size = 128;
    auto v = tomo::volume<D, T>(size);
    auto g = tomo::geometry::parallel<D, T>(v, size, size);
    auto f = tomo::modified_shepp_logan_phantom<T>(v);
    auto k = tomo::dim::joseph<D, T>(v);
    auto p = tomo::forward_projection<D, T>(f, g, k);

    auto x = tomo::reconstruction::sirt(v, g, k, p);
    tomo::ascii_plot(x);
}
```

### C++

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

### Python

The Python bindings expose the different concepts (images, volumes, geometries and dims) as well as the standard implemented algorithms.

# Building

## Dependencies

The following libraries are required:

External:
- [glm](http://glm.g-truc.net/0.9.8/index.html) header only mathematics library mimicking GLSL
- (optional) [MPI](http://www.nvidia.com/object/cuda_home_new.html) (>= 7.0)
- (optional) [RECAST3D](https://github.com/cicwi/RECAST3D) as an visualization server

Provided as submodules
- [Catch](https://github.com/philsquared/Catch), for unit tests
- [fmt](https://github.com/fmtlib/fmt) as an iostream replacement
- [bulk](https://github.com/jwbuurlage/Bulk) for distributed computing
- [cpptoml](https://github.com/skystrife/cpptoml) for reading specification and configuration files
- (optional) [zeromq](zeromq.org) for communicating with visualization servers
- (optional) [pybind11](https://github.com/pybind/pybind11) to generate Python bindings

The following build tools should be available on the system:
- [CMake](https://cmake.org/) (>= 3.0)
- Modern C++ compiler (with support for at least C++17), e.g. GCC >= 7.0 or clang >= 4.0

The library is being tested on Fedora 28, but is intended to be portable to other Linux distributions.

## Building process

The core of the library is header only, so it does not have to be built itself. We start with initializing the submodules:

```
git submodule init
git submodule update --remote
```

To build the examples:

```
cd build
cmake ..
make
```

The resulting binaries will be in the `bin` folder.

## Building the Python bindings

To generate the Python bindings:

```
cd python/build
cmake .
make
```

The Python bindings can be used through `tomo.py` in the `Python` folder, which also adds some rudimentary plotting functionality on top of the bindings.

## Building with optional features

To build the ZMQ and MPI based examples, run the following instead of `cmake .`.

```
cmake -DDISTRIBUTED=on -DPYTHON_BINDINGS=on .
```
