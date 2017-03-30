#include <string>

#include "fmt/format.h"

#include "tomo.hpp"
#include "util/simple_args.hpp"
#include "util/trees.hpp"

using T = float;
constexpr tomo::dimension D = 3_D;

int main(int argc, char* argv[]) {
    auto opts = tomo::options{argc, argv};

    // the usage should be:
    // ./partition -i geometry_meta_file -o output_partitioning -p procs -e epsilon
    if (!tomo::required_arguments(opts, {"-i", "-o", "-p", "-e"})) {
        std::cout << "Not all arguments passed.\n";
        return -1;
    }

    auto meta_file = tomo::arg(opts, "-i");
    auto output_file = tomo::arg(opts, "-o");
    auto processors = tomo::arg_as<int>(opts, "-p");
    auto epsilon = tomo::arg_as<T>(opts, "-e");

    auto problem = tomo::read_configuration<3_D, T>(meta_file);

    // now partition using problem.acquisition_geometry and
    // problem.object_volume
    auto tree = tomo::distributed::partition_bisection<D, T>(
        *problem.acquisition_geometry, problem.object_volume, processors,
        epsilon);
    auto neutral = tomo::to_neutral_tree<T>(tree, problem.object_volume);
    tomo::print_neutral_tree(neutral);

    // then save to:
    tomo::serialize_tree(neutral, output_file);

    return 0;
}
