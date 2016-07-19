#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <tomo.hpp>

int add(int i, int j) { return i + j; }

class image_2d_wrapper {
  public:
    image_2d_wrapper(int x, int y) : image_(tomo::volume<2_D>({x, y})) {}
    image_2d_wrapper(tomo::image<2_D> image) : image_(image) {}

    std::vector<double>& data() { return image_.data(); }

  private:
    tomo::image<2_D> image_;
};

image_2d_wrapper modified_shepp_logan_phantom(int x, int y) {
    return image_2d_wrapper(
        tomo::modified_shepp_logan_phantom<double>(tomo::volume<2_D>({x, y})));
}

void init_image(py::module& m) {
    m.def("add", &add, "A function which adds two numbers");

    py::class_<image_2d_wrapper>(m, "image")
        .def(py::init<int, int>())
        .def("data", &image_2d_wrapper::data,
             "return the underlying image data");

    m.def("modified_sl_phantom", &modified_shepp_logan_phantom, "obtain phantom");
}
