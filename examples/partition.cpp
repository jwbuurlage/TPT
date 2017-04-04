#include <string>

#include "fmt/format.h"

#include "tomo.hpp"
#include "util/plotter.hpp"
#include "util/simple_args.hpp"
#include "util/trees.hpp"

using T = float;
constexpr tomo::dimension D = 3_D;

int main(int argc, char* argv[]) {
    auto opts = tomo::options{argc, argv};

    // the usage should be:
    // ./partition -i geometry_meta_file -o output_partitioning -p procs -e
    // epsilon
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

    if (tomo::passed(opts, "--preview")) {
        tomo::util::ext_plotter<D, T> plotter(
            "tcp://localhost:5555", "Partitioning preview: " + output_file);
        auto part_bisected = bulk::tree_partitioning<D>(
            tomo::math::vec_to_array<D, int>(problem.object_volume.voxels()),
            processors, std::move(tree));
        plotter.send_partition_information(part_bisected, processors,
                                           problem.object_volume);

        auto proj_stack =
            tomo::projections<3_D, T>(*problem.acquisition_geometry);
        plotter.send_projection_data(
            (tomo::geometry::trajectory<3_D, T>&)(*problem
                                                       .acquisition_geometry),
            proj_stack, problem.object_volume);
    }

    return 0;
}
