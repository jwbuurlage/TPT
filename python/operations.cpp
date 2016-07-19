#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "tomo.hpp"

void init_operations(py::module& m) {
    m.def(
        "fb",
        &tomo::forward_projection<2_D, double, tomo::parallel_geometry<2_D, double>>,
        "forward project phantom");
}
