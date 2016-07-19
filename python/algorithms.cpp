#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "tomo.hpp"

using Geometry = tomo::parallel_geometry<2_D, double>;
using Projector = tomo::linear_projector<2_D, double>;

tomo::image<2_D, double>
art_wrapper(const tomo::volume<2_D>& v, const Geometry& g,
            const tomo::sinogram<2_D, double, Geometry, Projector>& p,
            double beta = 0.5) {
    return tomo::art(v, g, p, beta);
}

void init_algorithm(py::module& m) {
    m.def("art", &art_wrapper, "ART algorithm");
}
