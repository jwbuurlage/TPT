#include <cmath>
#include <iostream>

#include "tomo.hpp"

int main() {
    int k = 128;

    // create a 2D volume of size k x k
    auto v = tomo::volume<2_D>(k, k);

    auto f = tomo::modified_shepp_logan_phantom<double>(v);
    tomo::ascii_plot(f);

    // create a parallel geometry for the volume with k / 2 detectors and
    // k / 2 angles, for two dimensions
    auto g = tomo::parallel_geometry<2_D>(k / 2, k / 2, v);

    // simulate the experiment
    auto sino = tomo::forward_projection<2_D, double>(f, g);
    tomo::ascii_plot(sino);

    // run an algorithm to reconstruct the image
    auto x = tomo::art(v, g, sino);
    std::cout << "ART: " << std::endl;
    tomo::ascii_plot(x);

    // run an algorithm to reconstruct the image
    auto y = tomo::sart(v, g, sino);
    std::cout << "SART: " << std::endl;
    tomo::ascii_plot(y);

    // run an algorithm to reconstruct the image
    auto z = tomo::sirt(v, g, sino);
    std::cout << "SIRT: " << std::endl;
    tomo::ascii_plot(z);

    return 0;
}
