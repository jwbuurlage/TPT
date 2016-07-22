#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "tomo.hpp"

void init_operations(py::module& m) {
    py::class_<tomo::linear_projector<2_D, double>>(m, "linear_projector")
        .def(py::init<tomo::volume<2_D>>());

    m.def("forward_project",
          &tomo::forward_projection<2_D, double,
                                    tomo::parallel_geometry<2_D, double>>,
          "forward project phantom");
}
