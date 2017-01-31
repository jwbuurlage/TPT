/**
 * TODO:
 * - [ ] Arguments
 */
#include <cmath>
#include <iostream>

#include <glm/gtx/string_cast.hpp>

#include "tomo.hpp"
#include "util/plotter.hpp"

using T = double;

int main(int argc, char* argv[]) {
    // request a plotter scene
    auto plotter =
        tomo::ext_plotter<2_D, T>("tcp://localhost:5555", "Sequential test");

    auto opt = tomo::util::args(argc, argv);

    auto v = tomo::volume<2_D>(opt.k);
    auto f = tomo::modified_shepp_logan_phantom<T>(v);
    auto g = tomo::geometry::parallel<2_D, T>(opt.k, opt.k, v);

    auto kernel = tomo::dim::linear<2_D, T>(v);

    tomo::image<2_D, T> matrix(v);
    auto line = g.get_line(g.lines() / 2 + opt.k / 4);
    int k = 0;
    for (auto elem : kernel(line)) {
        matrix[elem.index] += elem.value;
        std::cout << k++ << " " << elem.index << ", " << elem.value << "\n";
    }

    plotter.plot(matrix);

    if (opt.sirt) {
       auto sino = tomo::forward_projection<2_D, T>(f, g, kernel);
       plotter.plot(sino.as_image());
        auto y = tomo::reconstruction::sirt(
            v, g, kernel, sino, 0.5, 10,
            {[&](tomo::image<2_D, T>& image) { plotter.plot(image); }});
    }
    return 0;
}
