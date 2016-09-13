#include <cmath>
#include <iostream>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include "fmt/format.h"

#include "tomo.hpp"

#include "../cuda/tomo_cuda.hpp"

int main(int argc, char* argv[]) {
	using T = float;

    int k = 0;
    int iterations = 0;
    T beta = 0.0f;

    po::options_description desc("Allowed arguments");
    desc.add_options()("help,h", "show help message")(
        "size,s", po::value<int>(&k)->default_value(256),
        "size of the phantom")("iterations,i",
                               po::value<int>(&iterations)->default_value(10),
                               "number of iterations")(
        "beta,b", po::value<T>(&beta)->default_value(0.5f),
        "value for update relaxation");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return -1;
    }

    // create a 2D volume of size k x k
    auto v = tomo::volume<2_D>(k, k);

    auto f = tomo::modified_shepp_logan_phantom<T>(v);
    tomo::ascii_plot(f);

    // create a parallel geometry for the volume with 250 detectors and
    // 180 angles, for two dimensions
    auto g = tomo::parallel_geometry<2_D, T>(k, k, v);

    // simulate the experiment
    auto sino = tomo::cuda::forward_projection<2_D, T>(f, g);
    // tomo::ascii_plot(sino);

    // run an algorithm to reconstruct the image
    auto x = tomo::cuda::sirt(v, g, sino, beta, iterations);

    fmt::print("CUDA SIRT:\n");
    tomo::ascii_plot(x);

    fmt::print("Parameters: size = {}x{}, iterations = {}, beta = {}\n", k, k,
               iterations, beta);

    return 0;
}
