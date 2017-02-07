#include <cmath>
#include <iostream>

#include "fmt/format.h"

#include "tomo.hpp"

using T = float;

template <tomo::dimension D>
void run(tomo::util::args opt) {
    // create a 2D volume of size k x k
    auto geom_and_volume =
        tomo::read_configuration<D, T>("data/geometries/parallel.toml");
    auto g = std::move(geom_and_volume.first);
    auto v = geom_and_volume.second;
    auto f = tomo::modified_shepp_logan_phantom<T>(v);

    tomo::ascii_plot(f);

    // simulate the experiment
    auto proj = tomo::dim::closest<D, T>(v);
    // auto proj = tomo::dim::joseph<T>(v);
    // auto proj = tomo::dim::closest<D, T>(v);
    auto sino = tomo::forward_projection<D, T>(f, *g, proj);
    // ascii_plot(sino);

    // run an algorithm to reconstruct the image
    if (opt.art) {
        auto x = tomo::reconstruction::art(v, *g, proj, sino, opt.beta,
                                           opt.iterations);
        fmt::print("ART\n");
        tomo::ascii_plot(x);
    }

    // run an algorithm to reconstruct the image
    if (opt.sart) {
        auto y = tomo::reconstruction::sart(v, *g, proj, sino, opt.beta,
                                            opt.iterations);
        fmt::print("SART\n");
        tomo::ascii_plot(y);
    }

    if (opt.sirt) {
        // run an algorithm to reconstruct the image
        auto z = tomo::reconstruction::sirt(v, *g, proj, sino, opt.beta,
                                            opt.iterations);
        fmt::print("SIRT\n");
        tomo::ascii_plot(z);
    }

    fmt::print("Parameters: size = {}x{}, iterations = {}, beta = {}\n",
               v.voxels()[0], v.voxels()[1], opt.iterations, opt.beta);
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
