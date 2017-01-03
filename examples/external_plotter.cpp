/**
 * TODO:
 * - [ ] Arguments
 */
#include <cmath>
#include <iostream>

#include "tomo.hpp"
#include "util/plotter.hpp"

using T = float;

int main() {
    // request a plotter scene
    auto plotter =
        tomo::ext_plotter<3_D, T>("tcp://localhost:5555", "Test plot");

    //auto args = tomo::util::args(argc, argv);

    int k = 92;

    auto v = tomo::volume<3_D>(k);
    auto f = tomo::modified_shepp_logan_phantom<T>(v);
    auto g = tomo::geometry::parallel<3_D, T>(k, k, v);

    auto proj = tomo::dim::closest<3_D, T>(v);
    auto sino = tomo::forward_projection<3_D, T>(f, g, proj);

    auto y = tomo::reconstruction::sirt(
        v, g, sino, 0.5, 10,
        {[&](tomo::image<3_D, T>& image) { plotter.plot(image); }});

    return 0;
}
