#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "tomo.hpp"

using Geometry = tomo::parallel_geometry<2_D, double>;
using Projector = tomo::linear_projector<2_D, double>;

tomo::image<2_D, double>
art_wrapper(const tomo::volume<2_D>& v, const Geometry& g,
            const tomo::sinogram<2_D, double, Geometry, Projector>& p,
            double beta = 0.5, int iterations = 10) {
    return tomo::art(v, g, p, beta, iterations);
}

tomo::image<2_D, double>
sart_wrapper(const tomo::volume<2_D>& v, const Geometry& g,
             const tomo::sinogram<2_D, double, Geometry, Projector>& p,
             double beta = 0.5, int iterations = 10) {
    return tomo::sart(v, g, p, beta, iterations);
}

tomo::image<2_D, double>
sirt_wrapper(const tomo::volume<2_D>& v, const Geometry& g,
             const tomo::sinogram<2_D, double, Geometry, Projector>& p,
             double beta = 0.5, int iterations = 10) {
    return tomo::sirt(v, g, p, beta, iterations);
}

tomo::image<2_D, double>
art_wrapper_initial(const tomo::volume<2_D>& v, const Geometry& g,
                    const tomo::sinogram<2_D, double, Geometry, Projector>& p,
                    double beta, tomo::image<2_D, double> initial) {
    return tomo::art(v, g, p, beta, 1,
                      optional<decltype(initial)>(initial));
}

tomo::image<2_D, double>
sart_wrapper_initial(const tomo::volume<2_D>& v, const Geometry& g,
                     const tomo::sinogram<2_D, double, Geometry, Projector>& p,
                     double beta, tomo::image<2_D, double> initial) {
    return tomo::sart(v, g, p, beta, 1,
                      optional<decltype(initial)>(initial));
}

tomo::image<2_D, double>
sirt_wrapper_initial(const tomo::volume<2_D>& v, const Geometry& g,
                     const tomo::sinogram<2_D, double, Geometry, Projector>& p,
                     double beta, tomo::image<2_D, double> initial) {
    return tomo::sirt(v, g, p, beta, 1,
                      optional<decltype(initial)>(initial));
}

void init_algorithm(py::module& m) {
    m.def("art", &art_wrapper, "ART reconstruction algorithm",
          py::arg("volume"), py::arg("geometry"), py::arg("projection"),
          py::arg("beta") = 0.5, py::arg("iterations") = 10);
    m.def("sart", &sart_wrapper, "SART reconstruction algorithm",
          py::arg("volume"), py::arg("geometry"), py::arg("projection"),
          py::arg("beta") = 0.5, py::arg("iterations") = 10);
    m.def("sirt", &sirt_wrapper, "SIRT reconstruction algorithm",
          py::arg("volume"), py::arg("geometry"), py::arg("projection"),
          py::arg("beta") = 0.5, py::arg("iterations") = 10);

    m.def("art_iteration", &art_wrapper_initial, "ART reconstruction algorithm",
          py::arg("volume"), py::arg("geometry"), py::arg("projection"),
          py::arg("beta"), py::arg("initial"));
    m.def("sart_iteration", &sart_wrapper_initial, "SART reconstruction algorithm",
          py::arg("volume"), py::arg("geometry"), py::arg("projection"),
          py::arg("beta"), py::arg("initial"));
    m.def("sirt_iteration", &sirt_wrapper_initial, "SIRT reconstruction algorithm",
          py::arg("volume"), py::arg("geometry"), py::arg("projection"),
          py::arg("beta"), py::arg("initial"));
}
