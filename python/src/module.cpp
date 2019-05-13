#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "python_common.hpp"
#include "tpt/tpt.hpp"

using namespace tpt::python;

namespace td = tpt::dim;
namespace tg = tpt::geometry;

auto left(bulk::util::binary_tree<bulk::util::split>::node* node) {
    return node->left.get();
}
auto right(bulk::util::binary_tree<bulk::util::split>::node* node) {
    return node->right.get();
}
auto root(bulk::util::binary_tree<bulk::util::split>* tree) {
    return tree->root.get();
}

void init_image(py::module& m) {
    py::class_<tpt::volume<2_D, T>>(m, "volume")
        .def(py::init<int>())
        .def("voxels", &tpt::volume<2_D, T>::voxels,
             "obtain a the number of voxels for each axis");

    py::class_<tpt::image<2_D, T>>(m, "image")
        .def(py::init<tpt::volume<2_D, T>>())
        .def("data", &tpt::image<2_D, T>::mutable_data,
             "obtain the underlying image data")
        .def("shape", &tpt::image<2_D, T>::shape, "obtain the image shape")
        .def("volume", &tpt::image<2_D, T>::get_volume,
             "obtain the underlying volume");

    m.def("modified_sl_phantom",
          py::overload_cast<const tpt::volume<2_D, T>&>(
              &tpt::modified_shepp_logan_phantom<T>),
          "create phantom");
}

void init_geometry(py::module& m) {
    // FIXME expose all vec{2,3,4}{i, f} with buffer interface
    py::class_<tpt::math::vec<2_D, int>>(m, "vec2i")
        .def(py::init<int, int>())
        .def_readwrite("x", &tpt::math::vec<2_D, int>::x)
        .def_readwrite("y", &tpt::math::vec<2_D, int>::y);

    py::class_<tpt::math::vec<3_D, int>>(m, "vec3i")
        .def(py::init<int, int, int>())
        .def_readwrite("x", &tpt::math::vec<3_D, int>::x)
        .def_readwrite("y", &tpt::math::vec<3_D, int>::y)
        .def_readwrite("z", &tpt::math::vec<3_D, int>::z);

    py::class_<tpt::math::vec<2_D, T>>(m, "vec2f")
        .def(py::init<T, T>())
        .def_readwrite("x", &tpt::math::vec<2_D, T>::x)
        .def_readwrite("y", &tpt::math::vec<2_D, T>::y);

    py::class_<tpt::math::vec<3_D, T>>(m, "vec3f")
        .def(py::init<T, T, T>())
        .def_readwrite("x", &tpt::math::vec<3_D, T>::x)
        .def_readwrite("y", &tpt::math::vec<3_D, T>::y)
        .def_readwrite("z", &tpt::math::vec<3_D, T>::z);

    py::class_<tpt::math::ray<2_D, T>>(m, "ray_2d");

    py::class_<tpt::dim::base<2_D, T>>(m, "base_kernel");
    py::class_<tpt::geometry::base<2_D, T>>(m, "base_geometry");

    py::class_<tpt::dim::linear<2_D, T>, tpt::dim::base<2_D, T>>(m, "linear")
        .def(py::init<tpt::volume<2_D, T>>());
    py::class_<tpt::dim::joseph<2_D, T>, tpt::dim::base<2_D, T>>(m, "joseph")
        .def(py::init<tpt::volume<2_D, T>>());
    py::class_<tpt::dim::closest<2_D, T>, tpt::dim::base<2_D, T>>(m,
                                                                    "closest")
        .def(py::init<tpt::volume<2_D, T>>());

    py::class_<tpt::geometry::parallel<2_D, T>, tpt::geometry::base<2_D, T>>(
        m, "parallel")
        .def(py::init<tpt::volume<2_D, T>, int>());

    py::class_<tpt::projections<2_D, T>>(m, "projections")
        .def("data", &tpt::projections<2_D, T>::mutable_data)
        .def("shape", &tpt::projections<2_D, T>::shape);
}

void init_geometry_3d(py::module& m) {
    py::class_<tpt::math::ray<3_D, T>>(m, "ray_3d")
        .def_readwrite("source", &tpt::math::ray<3_D, T>::source)
        .def_readwrite("detector", &tpt::math::ray<3_D, T>::detector);

    py::class_<tpt::volume<3_D, T>>(m, "volume_3d")
        .def(py::init<tpt::math::vec<3_D, int>>())
        .def(py::init<tpt::math::vec<3_D, int>, tpt::math::vec<3_D, T>,
                      tpt::math::vec<3_D, T>>())
        .def("voxels", &tpt::volume<3_D, T>::voxels,
             "obtain the voxel counts");

    py::class_<tpt::dim::base<3_D, T>>(m, "base_kernel_3d");

    py::class_<tpt::dim::linear<3_D, T>, tpt::dim::base<3_D, T>>(m,
                                                                   "linear_3d");
    py::class_<tpt::dim::joseph<3_D, T>, tpt::dim::base<3_D, T>>(m,
                                                                   "joseph_3d");
    py::class_<tpt::dim::closest<3_D, T>, tpt::dim::base<3_D, T>>(
        m, "closest_3d");

    py::class_<tpt::geometry::base<3_D, T>>(m, "base_geometry_3d");

    py::class_<tpt::geometry::cone_beam<T>, tpt::geometry::base<3_D, T>>(
        m, "cone_beam")
        .def(py::init<tpt::volume<3_D, T>, int, tpt::math::vec<2_D, T>,
                      tpt::math::vec<2_D, int>, T, T>());

    py::class_<tpt::geometry::projection<3_D, T>>(m, "projection")
        .def(py::init<tpt::math::vec<3_D, T>, tpt::math::vec<3_D, T>,
                      tpt::math::vec<3_D - 1, T>,
                      std::array<tpt::math::vec<3_D, T>, 3_D - 1>,
                      tpt::math::vec<3_D - 1, int>, bool>());

    py::class_<tpt::geometry::custom<T>, tpt::geometry::base<3_D, T>>(
        m, "custom_geometry")
        .def(py::init<std::vector<tpt::geometry::projection<3_D, T>>>());
}

void init_operations(py::module& m) {
    m.def("forward_project", &tpt::forward_projection<2_D, T>);
    m.def("back_project", &tpt::back_projection<2_D, T>);
    m.def("forward_project_3d", &tpt::forward_projection<3_D, T>);
    m.def("back_project_3d", &tpt::back_projection<3_D, T>);
}

void init_algorithm(py::module& m) {
    namespace tr = tpt::reconstruction;

    py::class_<std::function<void(tpt::image<2_D, T>&, int)>>(
        m, "callback_function")
        .def(py::init<>());

    m.def("art", &tr::art<2_D, T>, "ART reconstruction algorithm",
          py::arg("volume"), py::arg("geometry"), py::arg("kernel"),
          py::arg("projections"), py::arg("beta") = 0.5,
          py::arg("iterations") = 10,
          py::arg("callback") =
              std::function<void(tpt::image<2_D, T>&, int)>{});
    m.def("sart", &tr::sart<2_D, T>, "SART reconstruction algorithm",
          py::arg("volume"), py::arg("geometry"), py::arg("kernel"),
          py::arg("projections"), py::arg("beta") = 0.5,
          py::arg("iterations") = 10,
          py::arg("callback") =
              std::function<void(tpt::image<2_D, T>&, int)>{});
    m.def("sirt", &tr::sirt<2_D, T>, "SIRT reconstruction algorithm",
          py::arg("volume"), py::arg("geometry"), py::arg("kernel"),
          py::arg("projections"), py::arg("beta") = 0.5,
          py::arg("iterations") = 10,
          py::arg("callback") =
              std::function<void(tpt::image<2_D, T>&, int)>{},
          py::arg("box_constraint") = false, py::arg("box_min") = -1,
          py::arg("box_max") = 1);
}

void init_partitioning(py::module& m) {
    namespace td = tpt::distributed;

    py::class_<bulk::util::split>(m, "binary_tree_split")
        .def_readwrite("d", &bulk::util::split::d)
        .def_readwrite("a", &bulk::util::split::a);

    py::class_<bulk::util::binary_tree<bulk::util::split>::node>(
        m, "binary_tree_node")
        .def("left", &left, py::return_value_policy::reference)
        .def("right", &right, py::return_value_policy::reference)
        .def_readwrite(
            "value", &bulk::util::binary_tree<bulk::util::split>::node::value);

    py::class_<bulk::util::binary_tree<bulk::util::split>>(m,
                                                           "bulk_binary_tree")
        .def("root", &root, py::return_value_policy::reference);

    m.def("partition_bisection",
          &tpt::distributed::partition_bisection<3_D, T>,
          "GRCB partitioning of a geometry", py::arg("geometry"),
          py::arg("object_volume"), py::arg("processors"),
          py::arg("max_epsilon") = 0.2);
}

PYBIND11_MODULE(py_tpt, m) {
    m.doc() = "Bindings for TPT";

    init_image(m);
    init_geometry(m);
    init_geometry_3d(m);
    init_operations(m);
    init_algorithm(m);
    init_partitioning(m);
}
