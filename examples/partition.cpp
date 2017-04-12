#include <string>
#include <utility>
#include <vector>

#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

#include "fmt/format.h"

#include "tomo.hpp"
#include "util/plotter.hpp"
#include "util/simple_args.hpp"
#include "util/trees.hpp"

using T = float;
constexpr tomo::dimension D = 3_D;

void partition(std::string meta_file, std::string output_file, int processors,
               T epsilon, bool preview, tomo::util::report& table) {
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

    auto volume = problem.object_volume;
    volume.set_voxels({256, 256, 256});
    // TODO FIXME we need to upscale the geometry not the number of voxels..
    // its easier to load in large geometry, and downscale for partitioning,
    // then run on large

    auto large_tree = tomo::from_neutral_tree<T>(neutral, volume);

    auto part_bisected = bulk::tree_partitioning<D>(
        tomo::math::vec_to_array<D, int>(volume.voxels()),
        processors, std::move(large_tree));
    auto part_trivial = tomo::distributed::partition_trivial(
        *problem.acquisition_geometry, volume, processors);

    if (preview) {
        tomo::util::ext_plotter<D, T> plotter("tcp://localhost:5555",
                                              "PP: " + output_file);

        plotter.send_partition_information(part_bisected, processors,
                                           volume);

        auto proj_stack =
            tomo::projections<3_D, T>(*problem.acquisition_geometry);
        plotter.send_projection_data((tomo::geometry::trajectory<3_D, T>&)(*(
                                         problem.acquisition_geometry)),
                                     proj_stack, problem.object_volume);
    }

    // Store result to table
    auto overlap_trivial = tomo::distributed::communication_volume<D, T>(
        *problem.acquisition_geometry, volume, part_trivial);
    auto overlap_bisected = tomo::distributed::communication_volume<D, T>(
        *problem.acquisition_geometry, volume, part_bisected);
    T imp = (T)0.0;
    if (overlap_trivial != 0)
        imp = (overlap_trivial - overlap_bisected) / (T)overlap_trivial;

    auto name = fs::path(meta_file).stem();
    table.add_row(name);
    table.add_result(name, "trivial", overlap_trivial);
    table.add_result(name, "binary", overlap_bisected);
    table.add_result(name, "improvement", fmt::format("{:.1f}%", 100 * imp));
}

int main(int argc, char* argv[]) {
    auto opts = tomo::options{argc, argv};

    std::vector<std::pair<std::string, std::string>> ins_and_outs;
    int processors = 16;
    T epsilon = 0.05;

    // the usage should be:
    // ./partition -i geometry_meta_file -o output_partitioning -p procs -e
    // epsilon
    if (!(opts.required_arguments({"-i", "-o"}) ||
          opts.required_arguments({"--in", "--out"}))) {
        std::cout << "Input and/or output files not given.\n";
        return -1;
    }

    if (opts.required_arguments({"-i", "-o"})) {
        ins_and_outs.push_back(std::make_pair(opts.arg("-i"), opts.arg("-o")));
    } else if (opts.required_arguments({"--in", "--out"})) {
        auto ins = opts.args("--in");
        auto out_dir = opts.arg("--out");
        auto outs = ins;
        std::transform(outs.begin(), outs.end(), outs.begin(),
                       [&](std::string in) {
                           auto stem = fs::path(in).stem();
                           return out_dir + stem.string() + ".bsp";
                       });

        ins_and_outs = tomo::zip(ins, outs);
    }

    if (opts.passed("-p")) {
        processors = opts.arg_as<int>("-p");
    }
    if (opts.passed("-e")) {
        epsilon = opts.arg_as<T>("-e");
    }

    auto table = tomo::util::report(
        "Detector overlaps for geometric partitioning", "geometry");
    table.add_column("trivial");
    table.add_column("binary");
    table.add_column("improvement");

    for (auto in_and_out : ins_and_outs) {
        partition(in_and_out.first, in_and_out.second, processors, epsilon,
                  opts.passed("--preview"), table);
    }

    table.print();

    return 0;
}
