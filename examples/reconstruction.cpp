#include <cmath>
#include <iostream>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include "fmt/format.h"

#include "tomo.hpp"

int main(int argc, char* argv[]) {
    using T = float;

    int k = 0;
    int iterations = 0;
    float beta = 0.0f;

    // TODO: can we merge --use-cuda into these options, so that we dont need to
    // duplicate this
    // TODO: can we make option for which projector to use (now compile time,
    // need different options)
    po::options_description desc("Allowed arguments");
    desc.add_options()("help,h", "show help message")(
        "size,s", po::value<int>(&k)->default_value(256),
        "size of the phantom")("iterations,i",
                               po::value<int>(&iterations)->default_value(10),
                               "number of iterations")(
        "beta,b", po::value<float>(&beta)->default_value(0.5f),
        "value for update relaxation")("art", "reconstruct using ART")(
        "sart", "reconstruct using SART")("sirt", "reconstruct using SIRT");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return -1;
    }

    // TODO: want micro benchmarking merged here, instead of using `time ./prog`

    // create a 2D volume of size k x k
    auto v = tomo::volume<2_D>(k, k);

    auto f = tomo::modified_shepp_logan_phantom<T>(v);
    tomo::ascii_plot(f);

    // create a parallel geometry for the volume with 250 detectors and
    // 180 angles, for two dimensions
    auto g = tomo::parallel_geometry<2_D, T>(180, 250, v);

    // simulate the experiment
    // auto proj = tomo::linear_projector<2_D, T>(v);
    // auto proj = tomo::joseph_projector<T>(v);
    auto proj = tomo::closest_projector<T>(v);
    auto sino = tomo::forward_projection<2_D, T>(f, g, proj);

    // run an algorithm to reconstruct the image
    if (vm.count("art")) {
        auto x = tomo::art(v, g, sino, beta, iterations);
        fmt::print("ART: \n");
        tomo::ascii_plot(x);
    }

    // run an algorithm to reconstruct the image
    if (vm.count("sart")) {
        auto y = tomo::sart(v, g, sino, beta, iterations);
        fmt::print("SART: \n");
        tomo::ascii_plot(y);
    }

    if (vm.count("sirt")) {
        // run an algorithm to reconstruct the image
        auto z = tomo::sirt(v, g, sino, beta, iterations);
        fmt::print("SIRT: \n");
        tomo::ascii_plot(z);
    }

    fmt::print("Parameters: size = {}x{}, iterations = {}, beta = {}\n", k, k,
               iterations, beta);

    return 0;
}
