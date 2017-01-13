#pragma once

#include <boost/program_options.hpp>
namespace po = boost::program_options;

namespace tomo {
namespace util {

class args {
  public:
    int k = 32;
    int iterations = 10;
    float beta = 0.5f;

    bool two = true;
    bool three = false;

    bool art = false;
    bool sart = false;
    bool sirt = false;

    bool help = false;

    args(int argc, char* argv[]) { read_options_(argc, argv); }

  private:
    void read_options_(int argc, char* argv[]) {
        po::options_description desc("Allowed arguments");
        desc.add_options()("help,h", "show help message")(
            "size,s", po::value<int>(&k)->default_value(32),
            "size of the phantom")(
            "iterations,i", po::value<int>(&iterations)->default_value(10),
            "number of iterations")(
            "beta,b", po::value<float>(&beta)->default_value(0.5f),
            "value for update relaxation")("art", "reconstruct using ART")(
            "sart", "reconstruct using SART")("sirt", "reconstruct using SIRT")(
            "two", "do 2d reconstruction")("three", "do 3d reconstruction");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        art = vm.count("art");
        sart = vm.count("sart");
        sirt = vm.count("sirt");

        two = vm.count("two");
        three = vm.count("three");

        if (vm.count("help")) {
            std::cout << desc << "\n";
            help = true;
        }
    }
};

} // namespace util
} // namespace tomo
