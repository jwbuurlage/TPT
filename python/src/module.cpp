#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "python_common.hpp"
#include "tomo.hpp"

#ifdef USE_CUDA
#include "tomo_cuda.hpp"
#endif

using namespace tomo::python;

template <typename Sinogram>
std::array<int, 2> sino_dimensions(Sinogram& sino) {
    return std::array<int, 2>{sino.geometry().groups().x,
                              sino.geometry().groups().y};
}

template <typename G>
auto get_line(G& g, int i) {
    return g.get_line(i);
}

template <typename G>
auto lines(G& g) {
    return g.lines();
}

void init_image(py::module& m) {
    py::class_<tomo::volume<2_D>>(m, "volume").def(py::init<int, int>());

    py::class_<tomo::image<2_D, T>>(m, "image")
        .def(py::init<tomo::volume<2_D>>())
        .def("data", &tomo::image<2_D, T>::mutable_data,
             "obtain the underlying image data")
        .def("dimensions", &tomo::image<2_D, T>::dimensions,
             "obtain the image dimensions");

    m.def("modified_sl_phantom", &tomo::modified_shepp_logan_phantom<T>,
          "create phantom");
}

template <typename Ps, typename Gs>
void init_geometry(py::module& m, Ps ps, Gs gs) {
    py::class_<tomo::line<2_D, T>>(m, "line_2d");
    // TODO: add 'get_line' or some iterable adapter function
    hana::for_each(gs, [&](auto x) {
        using G = typename decltype(+x[1_c])::type;
        auto name = x[0_c] + "_geometry"s;
        using Init = typename decltype(
            hana::unpack(x[2_c], hana::template_<py::detail::init>))::type;
        py::class_<G>(m, name.c_str())
            .def(Init())
            .def("get_line", &get_line<G>, "obtain the i-th line");
    });

    auto combinations = hana::cartesian_product(hana::make_tuple(gs, ps));
    hana::for_each(combinations, [&](auto x) {
        using G = typename decltype(+(x[0_c][1_c]))::type;
        using P = typename decltype(+x[1_c][1_c])::type;
        using Sinogram = tomo::sinogram<2_D, T, G, P>;
        auto name = "sinogram_" + x[0_c][0_c] + "_"s + x[1_c][0_c];
        py::class_<Sinogram>(m, name.c_str())
            .def("data", &Sinogram::mutable_data,
                 "obtain the underlying image data")
            .def("dimensions", &sino_dimensions<Sinogram>,
                 "obtain the sinogram dimensions");
    });
}

template <typename Gs>
void init_geometry_3d(py::module& m, Gs gs) {
    py::class_<tomo::line<3_D, T>>(m, "line_3d")
        .def_readwrite("origin", &tomo::line<3_D, T>::origin)
        .def_readwrite("delta", &tomo::line<3_D, T>::delta)
        .def_readwrite("length", &tomo::line<3_D, T>::length);

    py::class_<tomo::math::vec<2_D, int>>(m, "vec2i").def(py::init<int, int>());
    py::class_<tomo::math::vec<3_D, T>>(m, "vec3f")
        .def(py::init<T, T, T>())
        .def_readwrite("x", &tomo::math::vec<3_D, T>::x)
        .def_readwrite("y", &tomo::math::vec<3_D, T>::y)
        .def_readwrite("z", &tomo::math::vec<3_D, T>::z);

    py::class_<tomo::volume<3_D>>(m, "volume_3d")
        .def(py::init<int, int, int>())
        .def("x", &tomo::volume<3_D>::x)
        .def("y", &tomo::volume<3_D>::y)
        .def("z", &tomo::volume<3_D>::z);

    hana::for_each(gs, [&](auto x) {
        using G = typename decltype(+x[1_c])::type;
        auto name = x[0_c] + "_geometry"s;
        using Init = typename decltype(
            hana::unpack(x[2_c], hana::template_<py::detail::init>))::type;
        py::class_<G>(m, name.c_str())
            .def(Init())
            .def("get_line", &get_line<G>, "obtain the i-th line")
            .def("lines", &lines<G>, "obtain the number of lines");
    });
}

template <typename Ps, typename Gs>
void init_operations(py::module& m, Ps ps, Gs gs) {
    hana::for_each(ps, [&](auto x) {
        auto name = x[0_c] + "_projector"s;
        using P = typename decltype(+x[1_c])::type;
        py::class_<P>(m, name.c_str()).def(py::init<tomo::volume<2_D>>());
    });

    auto combinations = hana::cartesian_product(hana::make_tuple(gs, ps));
    hana::for_each(combinations, [&](auto x) {
        using G = typename decltype(+(x[0_c][1_c]))::type;
        using P = typename decltype(+x[1_c][1_c])::type;
        m.def("forward_project", &tomo::forward_projection<2_D, T, G, P>);
    });
}

template <typename Ps, typename Gs>
void init_algorithm(py::module& m, Ps ps, Gs gs) {
    auto combinations = hana::cartesian_product(hana::make_tuple(gs, ps));
    hana::for_each(combinations, [&](auto x) {
        using G = typename decltype(+(x[0_c][1_c]))::type;
        using P = typename decltype(+x[1_c][1_c])::type;

        m.def("art", &tomo::reconstruction::art<2_D, T, G, P>,
              "ART reconstruction algorithm", py::arg("volume"),
              py::arg("geometry"), py::arg("projection"), py::arg("beta") = 0.5,
              py::arg("iterations") = 10);
        m.def("sart", &tomo::reconstruction::sart<2_D, T, G, P>,
              "SART reconstruction algorithm", py::arg("volume"),
              py::arg("geometry"), py::arg("projection"), py::arg("beta") = 0.5,
              py::arg("iterations") = 10);
        m.def("sirt", &tomo::reconstruction::sirt<2_D, T, G, P>,
              "SIRT reconstruction algorithm", py::arg("volume"),
              py::arg("geometry"), py::arg("projection"), py::arg("beta") = 0.5,
              py::arg("iterations") = 10);
    });
}

#ifdef USE_CUDA
template <typename Gs>
void init_cuda(py::module& m, Gs gs) {
    hana::for_each(gs, [&](auto x) {
        using G = typename decltype(+(x[1_c]))::type;
        m.def("cuda_forward_project",
              &tomo::cuda::forward_projection<2_D, T, G>);

        using Sinogram =
            tomo::sinogram<2_D, T, G, tomo::cuda::external_cuda_projector>;
        auto name = "sinogram_" + x[0_c] + "_cuda"s;
        py::class_<Sinogram>(m, name.c_str())
            .def("data", &Sinogram::mutable_data,
                 "obtain the underlying image data")
            .def("dimensions", &sino_dimensions<Sinogram>,
                 "obtain the sinogram dimensions");

        m.def("cuda_sart",
              &tomo::cuda::sart<2_D, T, G, tomo::cuda::external_cuda_projector>,
              "SART reconstruction algorithm using CUDA", py::arg("volume"),
              py::arg("geometry"), py::arg("projection"), py::arg("beta") = 0.5,
              py::arg("iterations") = 10);
    });
}
#endif

PYBIND11_PLUGIN(py_galactica) {
    py::module m("py_galactica", "bindings for galactica");

    // this is a list of the projector and geometry types, used to instantiate
    // the template algorithms and operations
    auto ps = hana::make_tuple(
        hana::make_tuple("linear"s, hana::type_c<tomo::dim::linear<2_D, T>>),
        hana::make_tuple("joseph"s, hana::type_c<tomo::dim::joseph<T>>),
        hana::make_tuple("closest"s, hana::type_c<tomo::dim::closest<2_D, T>>));

    // the third entry is the signature of the constructor
    auto gs = hana::make_tuple(
        hana::make_tuple("parallel"s,
                         hana::type_c<tomo::geometry::parallel<2_D, T>>,
                         hana::tuple_t<int, int, tomo::volume<2_D>>),
        hana::make_tuple("list"s, hana::type_c<tomo::geometry::list<2_D, T>>,
                         hana::tuple_t<std::vector<tomo::line<2_D, T>>>));

    auto gs3 = hana::make_tuple(
        hana::make_tuple("dual_axis_parallel"s,
                         hana::type_c<tomo::geometry::dual_axis_parallel<T>>,
                         hana::tuple_t<tomo::volume<3_D>, int, T,
                                       tomo::math::vec<2_D, int>>),
        hana::make_tuple("cone_beam"s,
                         hana::type_c<tomo::geometry::cone_beam<T>>,
                         hana::tuple_t<tomo::volume<3_D>, int, T,
                                       tomo::math::vec<2_D, int>, T, T>),
        hana::make_tuple("helical_cone_beam"s,
                         hana::type_c<tomo::geometry::helical_cone_beam<T>>,
                         hana::tuple_t<tomo::volume<3_D>, int, T,
                                       tomo::math::vec<2_D, int>, T, T>),
        hana::make_tuple("dynamic_cone_beam"s,
                         hana::type_c<tomo::geometry::dynamic_cone_beam<T>>,
                         hana::tuple_t<tomo::volume<3_D>, int, T,
                                       tomo::math::vec<2_D, int>>),
        hana::make_tuple("laminography"s,
                         hana::type_c<tomo::geometry::laminography<T>>,
                         hana::tuple_t<tomo::volume<3_D>, int, T,
                                       tomo::math::vec<2_D, int>, T, T, T, T>),
        hana::make_tuple("tomosynthesis"s,
                         hana::type_c<tomo::geometry::tomosynthesis<T>>,
                         hana::tuple_t<tomo::volume<3_D>, int, T,
                                       tomo::math::vec<2_D, int>, T, T, T>)
        );

    init_image(m);
    init_geometry(m, ps, gs);
    init_geometry_3d(m, gs3);
    init_operations(m, ps, gs);
    init_algorithm(m, ps, gs);

#ifdef USE_CUDA
    init_cuda(m, gs);
#endif

    return m.ptr();
}
