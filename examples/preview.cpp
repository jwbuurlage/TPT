#include <cmath>
#include <fstream>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

using namespace std::string_literals;

#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

#include "fmt/format.h"

#include "tomos/tomos.hpp"
#include "tomos/util/plotter.hpp"
#include "tomos/util/simple_args.hpp"
#include "tomos/util/trees.hpp"

using T = float;

void preview(std::string tree_file, std::string geom_file, int k, int p) {
    auto obj_vol = tomo::volume<3_D, T>(k);
    auto name = fs::path(tree_file).stem().string();
    auto tree_part = tomo::load_partitioning(tree_file, obj_vol, log2(p));

    auto problem = tomo::read_configuration<3_D, T>(geom_file, k);

    tomo::util::ext_plotter<3_D, T> plotter("tcp://localhost:5555", name);
    plotter.send_geometry(*problem.acquisition_geometry, obj_vol);
    plotter.send_partition_information(*tree_part, p, obj_vol);

    auto weights = tomo::distributed::voxel_weights<3_D, T>(*problem.acquisition_geometry, obj_vol);
    plotter.plot(weights);
}

int main(int argc, char* argv[]) {
    auto opts = tomo::options{argc, argv};

    
    if (!(opts.required_arguments({"-t", "-g", "-k", "-p"}))) {
        std::cout << "Required: -t TREE -g GEOM -k SIZE -p PROCS\n";
        return -1;
    }

    preview(opts.arg("-t"), opts.arg("-g"), opts.arg_as<int>("-k"), opts.arg_as<int>("-p"));

    return 0;
}
