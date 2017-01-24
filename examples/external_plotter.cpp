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

    auto proj = tomo::dim::joseph<2_D, T>(v);

    tomo::image<2_D, T> matrix(v);
    auto line = g.get_line(g.lines() / 2 + opt.k / 4);
    for (auto elem : proj(line)) {
        matrix[elem.index] += elem.value;
    }

    plotter.plot(matrix);

    if (opt.sirt) {
       auto sino = tomo::forward_projection<2_D, T>(f, g, proj);
       plotter.plot(sino.as_image());
        auto y = tomo::reconstruction::sirt(
            v, g, sino, 0.5, 10,
            {[&](tomo::image<2_D, T>& image) { plotter.plot(image); }});
    }
    return 0;
}
