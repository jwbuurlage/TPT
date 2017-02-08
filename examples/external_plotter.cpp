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
constexpr tomo::dimension D = 3_D;

int main(int argc, char* argv[]) {
    // request a plotter scene
    auto plotter =
        tomo::ext_plotter<D, T>("tcp://localhost:5555", "Sequential test");

    auto opt = tomo::util::args(argc, argv);

    auto geom_and_volume =
        tomo::read_configuration<D, T>("data/geometries/parallel.toml");
    auto geom = std::move(geom_and_volume.first);
    auto v = geom_and_volume.second;
    auto f = tomo::modified_shepp_logan_phantom<T>(v);

    auto& g = *geom;

    std::cout << glm::to_string(v.origin()) << "\n";
    std::cout << glm::to_string(v.voxels()) << "\n";

    auto kernel = tomo::dim::joseph<D, T>(v);

//    tomo::image<D, T> matrix(v);
//    auto line = g.get_line(g.lines() / 2 + opt.k / 4);
//    for (auto elem : kernel(line)) {
//        matrix[elem.index] += elem.value;
//    }

    plotter.plot(f);
    tomo::ascii_plot(f);

    if (opt.sirt) {
        auto sino = tomo::forward_projection<D, T>(f, g, kernel);
        auto y = tomo::reconstruction::sirt(
            v, g, kernel, sino, 0.5, 10,
            {[&](tomo::image<D, T>& image) { plotter.plot(image); }});
    }

    return 0;
}
