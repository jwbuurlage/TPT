#include <pybind11/pybind11.h>
namespace py = pybind11;

void init_image(py::module&);

PYBIND11_PLUGIN(tomography) {
    py::module m("tomography", "bindings for galactica");

    init_image(m);

    return m.ptr();
}
