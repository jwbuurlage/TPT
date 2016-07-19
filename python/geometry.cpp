#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "tomo.hpp"

using Sinogram =
    tomo::sinogram<2_D, double, tomo::parallel_geometry<2_D, double>,
                   tomo::linear_projector<2_D, double>>;

std::array<int, 2> sino_dimensions(Sinogram& sino) {
    return std::array<int, 2>{sino.geometry().groups().x,
                              sino.geometry().groups().y};
}

void init_geometry(py::module& m) {
    py::class_<tomo::parallel_geometry<2_D, double>>(m, "parallel_geometry")
        .def(py::init<int, int, tomo::volume<2_D>>());

    py::class_<Sinogram>(m, "sinogram")
        .def("data", &Sinogram::mutable_data,
             "obtain the underlying image data")
        .def("dimensions", &sino_dimensions, "obtain the sinogram dimensions");
}
