# Tomos

Tomos is a library with many tools for tomographic reconstruction. Its features over other 'tomography toolboxes' include:

- Written in modern C++, with OO / generic programming style like STL
- Tomos supports a general model for distributed computing, enabling the resulting algorithms to run on clusters
- It is not bound to CUDA, MATLAB, or Python, although Python bindings and GPU computing are supported
- The CPU implementations are efficient
- Support for 2D, 3D or higher dimensional reconstructions.

# Dependencies

We use:
- Catch
- fmtlib/fmt
- glm
- boost::program_options
- thrust

For the Python bindings:
- pybind11
- boost::hana

To build:
- CMake (>= 3.0)
- Modern compiler (with support for at least C++14)

For GPU support:
- CUDA (>= 7.0)
