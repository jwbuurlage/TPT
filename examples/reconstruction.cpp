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

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include "fmt/format.h"

#include "tomo.hpp"

using T = float;

struct options {
    int k = 0;
    int iterations = 0;
    float beta = 0.0f;

    bool two = true;
    bool three = false;

    bool art = false;
    bool sart = false;
    bool sirt = false;

    bool help = false;
};

options read_options(int argc, char* argv[]) {
    options opt;

    po::options_description desc("Allowed arguments");
    desc.add_options()("help,h", "show help message")(
        "size,s", po::value<int>(&opt.k)->default_value(32),
        "size of the phantom")(
        "iterations,i", po::value<int>(&opt.iterations)->default_value(10),
        "number of iterations")(
        "beta,b", po::value<float>(&opt.beta)->default_value(0.5f),
        "value for update relaxation")("art", "reconstruct using ART")(
        "sart", "reconstruct using SART")("sirt", "reconstruct using SIRT")(
        "two", "do 2d reconstruction")("three", "do 3d reconstruction");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    opt.art = vm.count("art");
    opt.sart = vm.count("sart");
    opt.sirt = vm.count("sirt");

    opt.two = vm.count("two");
    opt.three = vm.count("three");

    if (vm.count("help")) {
        std::cout << desc << "\n";
        opt.help = true;
    }

    // TODO: can we merge --use-cuda into these options, so that we dont need to
    // duplicate this

    return opt;
}

template <tomo::dimension D>
void run(options opt) {
    // create a 2D volume of size k x k
    auto v = tomo::volume<D>(opt.k);

    auto f = tomo::modified_shepp_logan_phantom<T>(v);

    tomo::ascii_plot(f);

    // auto g = tomo::geometry::helical_cone_beam<T>(v, k, (T)1.0, {k, k});
    auto g = tomo::geometry::parallel<D, T>(opt.k, opt.k, v);

    //    // simulate the experiment
    auto proj = tomo::dim::closest<D, T>(v);
    //    // auto proj = tomo::dim::joseph<T>(v);
    //    auto proj = tomo::dim::closest<D, T>(v);
    auto sino = tomo::forward_projection<D, T>(f, g, proj);
    //ascii_plot(sino);

    // run an algorithm to reconstruct the image
    if (opt.art) {
        auto x =
            tomo::reconstruction::art(v, g, sino, opt.beta, opt.iterations);
        fmt::print("ART: \n");
        tomo::ascii_plot(x);
    }

    // run an algorithm to reconstruct the image
    if (opt.sart) {
        auto y =
            tomo::reconstruction::sart(v, g, sino, opt.beta, opt.iterations);
        fmt::print("SART: \n");
        tomo::ascii_plot(y);
    }

    if (opt.sirt) {
        // run an algorithm to reconstruct the image
        auto z =
            tomo::reconstruction::sirt(v, g, sino, opt.beta, opt.iterations);
        fmt::print("SIRT: \n");
        tomo::ascii_plot(z);
    }

    fmt::print("Parameters: size = {}x{}, iterations = {}, beta = {}\n", opt.k,
               opt.k, opt.iterations, opt.beta);
}

int main(int argc, char* argv[]) {
    auto opt = read_options(argc, argv);

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
