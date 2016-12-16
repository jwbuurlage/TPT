/**
 * TODO:
 * - [ ] 
 */
#include <cmath>
#include <iostream>

#include "tomo.hpp"
#include "util/plotter.hpp"

using T = float;

int main() {
    // request a plotter scene
    auto plotter =
        tomo::ext_plotter<2_D, T>("tcp://localhost:5555", "Test plot");

    int k = 512;

    auto v = tomo::volume<2_D>(k);
    auto f = tomo::modified_shepp_logan_phantom<T>(v);
    auto g = tomo::geometry::parallel<2_D, T>(k, k, v);

    auto proj = tomo::dim::closest<2_D, T>(v);
    auto sino = tomo::forward_projection<2_D, T>(f, g, proj);

    auto y = tomo::reconstruction::sirt(
        v, g, sino, 0.5, 10,
        {[&](tomo::image<2_D, T>& image) { plotter.plot(image); }});

    return 0;
}
