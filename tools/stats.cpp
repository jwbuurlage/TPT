#include <filesystem>
namespace fst = std::filesystem;

#include "bulk/bulk.hpp"
#include "tpt/tpt.hpp"
#include <CLI/CLI.hpp>

using T = float;
namespace td = tpt::distributed;

int main(int argc, char** argv) {
    CLI::App app{"Partitioning stats for distributed tomography"};

    int p = 2;
    int k = p * 8;
    std::vector<std::string> geoms;
    std::string tree_dir = ".";

    // basic options for partitioning
    app.add_option("-p", p, "number of processors to partition for", p);
    app.add_option("-k", k, "number of voxels for stats", k);

    // file input / output
    app.add_option("--geometries", geoms,
                   "the input geometries (.toml files)")
      ->required();
    app.add_option("--output-dir", tree_dir,
                   "the output directory for partitionings (.bsp files)",
                   tree_dir.c_str());

    CLI11_PARSE(app, argc, argv);

    for (auto geom : geoms) {
        // get base name and partitioning file
        auto name = fst::path(geom).stem().string();
        auto tree_file = tree_dir + "/" + name + ".bsp";

        // load problem
        auto problem = tpt::read_configuration<3_D, T>(geom, k);
        auto& g = *problem.acquisition_geometry;
        auto v = problem.object_volume;

        // load partitioning from file
        auto tree_part = tpt::load_partitioning(tree_file, v, log2(p));

        // compute stats
        auto imbalance = td::load_imbalance(v, *tree_part, g);
        auto comvol = td::communication_volume<3_D, T>(g, v, *tree_part);

        // compute stats
        auto regcomvol = td::regularizer_volume(v, *tree_part);

        // output stats
        std::cout << name << ": " << comvol << " + " << regcomvol << " = "
                  << comvol + regcomvol << " (" << imbalance << ")\n";
    }
}
