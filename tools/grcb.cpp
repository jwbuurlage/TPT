#include <CLI/CLI.hpp>

#include <filesystem>
#include <string>
namespace fst = std::filesystem;
using namespace std::string_literals;

using T = float;

#include "bulk/bulk.hpp"

#include "tpt/tpt.hpp"
namespace td = tpt::distributed;

int main(int argc, char** argv) {
    CLI::App app{"GRCB partitioner for distributed tomography"};

    int p = 2;
    int k = p * 8;
    std::vector<std::string> geometries;
    std::string partitioning_dir = ".";

    // basic options for partitioning
    app.add_option("-p", p, "number of processors to partition for", p);
    app.add_option("-k", k, "number of voxels for stats", k);

    // file input / output
    app.add_option("--geometries", geometries,
                   "the input geometries (.toml files)")
        ->required();
    app.add_option("--output-dir", partitioning_dir,
                   "the output directory for partitionings (.bsp files)",
                   partitioning_dir.c_str());

    // set partitioning strategy
    auto strategies = std::map<std::string, tpt::grcb::flag>{
        {"midway", tpt::grcb::flag::split_midway},
        {"smart", tpt::grcb::flag::split_smart},
        {"simple", tpt::grcb::flag::split_smart_simple},
        {"load", tpt::grcb::flag::split_equal_load},
    };
    tpt::grcb::flag strategy = tpt::grcb::flag::split_smart_simple;
    auto strategy_opt =
        app.add_option("--strategy", strategy, "load balancing strategy to use",
                       "simple")
            ->transform(CLI::CheckedTransformer(strategies));

    // set load balancing precision
    auto precisions = std::map<std::string, tpt::grcb::flag>{
        {"low", tpt::grcb::flag::precision_low},
        {"medium", tpt::grcb::flag::precision_medium},
        {"high", tpt::grcb::flag::precision_high},
    };
    tpt::grcb::flag level = tpt::grcb::flag::precision_medium;
    app.add_option("--precision", level, "load balancing accuracy", "medium")
        ->transform(CLI::CheckedTransformer(precisions));

    // set integration strategy
    auto kinds = std::map<std::string, tpt::grcb::flag>{
        {"2d", tpt::grcb::flag::integration_two_dim},
        {"3d", tpt::grcb::flag::integration_three_dim},
        {"cone", tpt::grcb::flag::integration_cone},
    };
    auto kind = tpt::grcb::flag::integration_cone;
    app.add_option("--integration", kind,
                   "integration method for 'load' strategy", "cone")
        ->transform(CLI::CheckedTransformer(kinds))
        ->needs(strategy_opt);

    // parse the CLI
    CLI11_PARSE(app, argc, argv);

    // if partitioning output dir does not exist, make it
    fst::create_directories(partitioning_dir);

    auto flags = std::vector<tpt::grcb::flag>{kind, level, strategy};

    // next we actually run the experiment

    auto report = bulk::util::table("Results"s, "geometry");
    report.columns("V", "epsilon", "part time", "eval time");

    for (auto geometry_file : geometries) {
        auto name = fst::path(geometry_file).stem().string();
        auto problem = tpt::read_configuration<3_D, T>(geometry_file, k);

        auto v = problem.object_volume;
        auto& g = *problem.acquisition_geometry;
        auto dt = bulk::util::timer();
        auto root = tpt::grcb::partition(v, g, std::log2(p), flags);
        auto part_time = dt.get();

        auto outfile = fst::path(partitioning_dir) / name;
        outfile += ".bsp";

        {
            auto of = std::ofstream(outfile);
            tpt::grcb::print_tree(root, of);

            std::cout << "Wrote: " << outfile << " to disk.\n";
        }

        auto dt2 = bulk::util::timer();
        auto tree_part = tpt::load_partitioning(outfile, v, log2(p));
        auto imbalance = td::load_imbalance(v, *tree_part, g);
        auto comvol = td::communication_volume<3_D, T>(g, v, *tree_part);
        auto eval_time = dt2.get();
        report.row(name, comvol, imbalance, part_time, eval_time);
    }

    std::cout << report.print() << "\n";

    return 0;
}
