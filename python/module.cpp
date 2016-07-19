#include <pybind11/pybind11.h>
namespace py = pybind11;

void init_image(py::module&);
void init_geometry(py::module&);
void init_operations(py::module&);
void init_algorithm(py::module&);

PYBIND11_PLUGIN(py_galactica) {
    py::module m("py_galactica", "bindings for galactica");

    init_image(m);
    init_geometry(m);
    init_operations(m);
    init_algorithm(m);

    return m.ptr();
}
