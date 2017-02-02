/**
 * TODO:
 * - [x] merge 2d and 3d reconstruction example into a single file
 * - [ ] allow choosing geometry for reconstruction
 * - [ ] fix sinogram for 3d (stack of projections != sinogram)
 * - [ ] fix plotting of sinogram
 * - [ ] see where reconstruction algorithms fail or get slow
 * - [ ] read in tiff stack and apply to real data
 * - [ ] optimize everything
 */
#include <cmath>
#include <iostream>

#include "fmt/format.h"

#include "tomo.hpp"

using T = float;

template <tomo::dimension D>
void run(tomo::util::args opt) {
    // create a 2D volume of size k x k
    auto v = tomo::volume<D, T>(opt.k);

    auto f = tomo::modified_shepp_logan_phantom<T>(v);

    tomo::ascii_plot(f);

    // auto g = tomo::geometry::helical_cone_beam<T>(v, k, (T)1.0, {k, k});
    auto g = tomo::geometry::parallel<D, T>(v, opt.k, opt.k);

    // simulate the experiment
    auto proj = tomo::dim::closest<D, T>(v);
    // auto proj = tomo::dim::joseph<T>(v);
    // auto proj = tomo::dim::closest<D, T>(v);
    auto sino = tomo::forward_projection<D, T>(f, g, proj);
    // ascii_plot(sino);

    // run an algorithm to reconstruct the image
    if (opt.art) {
        auto x = tomo::reconstruction::art(v, g, proj, sino, opt.beta,
                                           opt.iterations);
        fmt::print("ART\n");
        tomo::ascii_plot(x);
    }

    // run an algorithm to reconstruct the image
    if (opt.sart) {
        auto y = tomo::reconstruction::sart(v, g, proj, sino, opt.beta,
                                            opt.iterations);
        fmt::print("SART\n");
        tomo::ascii_plot(y);
    }

    if (opt.sirt) {
        // run an algorithm to reconstruct the image
        auto z = tomo::reconstruction::sirt(v, g, proj, sino, opt.beta,
                                            opt.iterations);
        fmt::print("SIRT\n");
        tomo::ascii_plot(z);
    }

    fmt::print("Parameters: size = {}x{}, iterations = {}, beta = {}\n", opt.k,
               opt.k, opt.iterations, opt.beta);
}

int main(int argc, char* argv[]) {
    auto opt = tomo::util::args(argc, argv);

    // TODO: can we make option for which projector to use (now compile time,
    // need different options)
    // TODO: want micro benchmarking merged here, instead of using `time ./prog`

    if (opt.two) {
        run<2_D>(opt);
    }
    if (opt.three) {
        run<3_D>(opt);
    }

    return 0;
}
