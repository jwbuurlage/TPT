#include <cmath>
#include <fstream>
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

namespace td = tomo::distributed;

using T = float;

void compute(std::string tree_file, std::string geometry_file, int a, int b,
             tomo::util::report& table, bool preview) {
    (void)geometry_file;

    auto name = fs::path(tree_file).stem().string();
    table.add_row(name + " V_b"s);
    table.add_row(name + " V_t"s);
    table.add_row(name + " g"s);
    table.add_row(name + " eps"s);
    for (int p = a; p <= b; p *= 2) {
        int k = tomo::math::min(p * 4, 512);
        auto obj_vol = tomo::volume<3_D, T>({k, k, k}, {0, 0, 0}, {1, 1, 1});
        auto geom = tomo::geometry::cone_beam<T>(obj_vol, 1, {1, 1}, {1, 1},
                                                 1.0f, 1.0f);
        auto tree_part = tomo::load_partitioning(tree_file, obj_vol, log2(p));
        auto problem = tomo::read_configuration<3_D, T>(geometry_file, k);
        auto& geometry = *problem.acquisition_geometry;

        auto part_trivial = tomo::distributed::partition_trivial(
            *problem.acquisition_geometry, obj_vol, p);

        auto overlap_bisected =
            td::communication_volume<3_D, T>(geometry, obj_vol, *tree_part);
        auto overlap_trivial =
            td::communication_volume<3_D, T>(geometry, obj_vol, part_trivial);

        T imp = (T)0.0;
        if (overlap_trivial != 0)
            imp = (overlap_trivial - overlap_bisected) / (T)overlap_trivial;

        auto imbalance =
            tomo::distributed::load_imbalance(obj_vol, *tree_part, geometry);

        table.add_result(name + " V_b"s, "p = "s + std::to_string(p),
                         overlap_bisected);
        table.add_result(name + " V_t"s, "p = "s + std::to_string(p),
                         overlap_trivial);
        table.add_result(name + " g"s, "p = "s + std::to_string(p),
                         fmt::format("{:.1f}%", 100 * imp));
        table.add_result(name + " eps"s, "p = "s + std::to_string(p),
                         fmt::format("{:.2f}", imbalance));

        // auto partitioning = from_obj_and_tree(tree, obj_vol, p);
        // split in 'i', make partitioning, compute comvol, add gain and trivial
        // to table
        // this comes back to the whole restricting geometries.. or it takes
        // forever, maybe just first run a test with small geometries based on
        // partitionings by the big geometries
        // we want the detector shape to be maybe 2 * p
        // number of projections also 2 * p
        if (preview) {
            tomo::util::ext_plotter<3_D, T> plotter(
                "tcp://localhost:5555",
                "PP: "s + tree_file + " p = "s + std::to_string(p));
            plotter.send_partition_information(*tree_part, p, obj_vol);
            plotter.send_projection_data(
                (tomo::geometry::trajectory<3_D, T>&)(geometry),
                problem.projection_stack, obj_vol);
        }
    }
}

void usage(std::string program_name) {
    std::cout
        << "USAGE: " << program_name
        << " --in TREES --geom GEOM -o TABLE_FILE -a MIN_PROC -b MAX_PROC\n";
}

int main(int argc, char* argv[]) {
    auto opts = tomo::options{argc, argv};

    // the usage should be:
    // ./genvol --in [trees] --geom [geom] -o [file] -a [minproc] -b [maxproc]
    if (!(opts.required_arguments({"--in", "--geom", "-o", "-a", "-b"}))) {
        usage(argv[0]);
        return -1;
    }

    auto ins = opts.args("--in");
    auto geoms = opts.args("--geom");
    auto out = opts.arg("-o");
    auto a = opts.arg_as<int>("-a");
    auto b = opts.arg_as<int>("-b");

    auto table = tomo::util::report("Communication volumes", "geometry");
    for (int i = a; i <= b; i *= 2) {
        table.add_column("p = "s + std::to_string(i));
    }

    assert(ins.size() == geoms.size());

    for (auto i = 0u; i < ins.size(); ++i) {
        std::cout << ins[i] << " " << geoms[i] << "\n";
        compute(ins[i], geoms[i], a, b, table, opts.passed("--preview"));
    }

    table.print();

    std::ofstream of(out, std::ios::out);
    table.print(of);

    return 0;
}
