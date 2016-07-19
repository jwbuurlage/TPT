#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "tomo.hpp"

void init_image(py::module& m) {
    py::class_<tomo::volume<2_D>>(m, "volume").def(py::init<int, int>());

    py::class_<tomo::image<2_D, double>>(m, "image")
        .def(py::init<tomo::volume<2_D>>())
        .def("data", &tomo::image<2_D, double>::mutable_data,
             "obtain the underlying image data")
        .def("dimensions", &tomo::image<2_D, double>::dimensions,
             "obtain the image dimensions");

    m.def("modified_sl_phantom", &tomo::modified_shepp_logan_phantom<double>,
          "create phantom");
}
