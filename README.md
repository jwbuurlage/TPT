# Galactica

Personal goal: convince that ASTRA 2.0 is needed with the following features:

* Written in modern C++, with OO / generic programming style like STL
* Distributed from the start
* Not bound to CUDA / MATLAB / Python, should have no reference anywhere in the main code
* Better build system (CMake)
* Data oriented => faster
* No or very little specific code for 2D or 3D in main code

# Dependencies

We use:

* Catch
* fmtlib/fmt
* glm
* boost::hana
* boost::program_options
* pybind11
* thrust

Build

* cmake
* cuda
* modern compiler
