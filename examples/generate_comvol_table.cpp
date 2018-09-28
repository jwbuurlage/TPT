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

namespace td = tomo::distributed;

using T = float;

std::mutex g_result_mutex;

void compute(std::string tree_dir, std::string geometry_file, int a, int b,
             tomo::util::report& table, bool preview) {

    auto name = fs::path(geometry_file).stem().string();
    table.add_row(name + " V_b"s);
    table.add_row(name + " V_cube"s);
    table.add_row(name + " V_b (reg)"s);
    table.add_row(name + " V_t"s);
    table.add_row(name + " g"s);
    table.add_row(name + " eps"s);
    table.add_row(name + " eps_triv"s);
    for (int p = a; p <= b; p *= 2) {
        // use tree dir and geometry file to  get tree file
        // ...
        auto tree_file = tree_dir + "/" + name + ".bsp";

        int k = tomo::math::min(512, p * 4);
        // auto obj_vol = tomo::volume<3_D, T>({k, k, k}, {0, 0, 0}, {1, 1, 1});
        auto problem = tomo::read_configuration<3_D, T>(geometry_file, k);
        auto& geometry = *problem.acquisition_geometry;
        auto obj_vol = problem.object_volume;
        auto tree_part = tomo::load_partitioning(tree_file, obj_vol, log2(p));
        tomo::print_tree(tree_part->splits());

        auto overlap_bisected =
            td::communication_volume<3_D, T>(geometry, obj_vol, *tree_part);

        auto main_d = tree_part->splits().root->value.d;
        auto part_trivial = bulk::block_partitioning<3_D, 1>(
            tomo::math::vec_to_array<3_D, int>(obj_vol.voxels()), {p},
            {main_d});

        auto overlap_trivial =
            td::communication_volume<3_D, T>(geometry, obj_vol, part_trivial);

        T imp = (T)0.0;
        if (overlap_trivial != 0)
            imp = (overlap_trivial - overlap_bisected) / (T)overlap_trivial;

        auto imbalance =
            tomo::distributed::load_imbalance(obj_vol, *tree_part, geometry);
        auto triv_imbalance =
            tomo::distributed::load_imbalance(obj_vol, part_trivial, geometry);

        auto comvol_reg = td::regularizer_volume<3_D, T>(obj_vol, *tree_part);

        auto cube_vol = 0;

        if (p == 64) {
            auto part_cube = bulk::block_partitioning<3_D, 3_D>(
                tomo::math::vec_to_array<3_D, int>(obj_vol.voxels()),
                {4, 4, 4});
            cube_vol =
                td::communication_volume<3_D, T>(geometry, obj_vol, part_cube);
        }

        {
            std::lock_guard<std::mutex> guard(g_result_mutex);
            table.add_result(name + " V_b"s, "p = "s + std::to_string(p),
                             overlap_bisected);
            table.add_result(name + " V_t"s, "p = "s + std::to_string(p),
                             overlap_trivial);
            table.add_result(name + " V_b (reg)"s, "p = "s + std::to_string(p),
                             comvol_reg);

            if (p == 64) {
                table.add_result(name + " V_cube"s, "p = "s + std::to_string(p),
                                 cube_vol);
            }

            table.add_result(name + " g"s, "p = "s + std::to_string(p),
                             fmt::format("{:.1f}%", 100 * imp));
            table.add_result(name + " eps"s, "p = "s + std::to_string(p),
                             fmt::format("{:.2f}", imbalance));
            table.add_result(name + " eps_triv"s, "p = "s + std::to_string(p),
                             fmt::format("{:.2f}", triv_imbalance));

            std::cout << tree_file + " imbalance: " << imbalance << "\n";
        }

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
            //            plotter.send_projection_data(
            //                (tomo::geometry::trajectory<3_D, T>&)(geometry),
            //                problem.projection_stack, obj_vol);
        }
    }
}

void usage(std::string program_name) {
    std::cout << "USAGE: " << program_name
              << " --in TREE_MAIN_FOLDER --geom "
                 "GEOMS -o TABLE_FILE -a MIN_PROC "
                 "-b MAX_PROC\n";
}

int main(int argc, char* argv[]) {
    auto opts = tomo::options{argc, argv};

    // the usage should be:
    // ./comvol --in [trees] --geom [geom] -o [file] -a [minproc] -b [maxproc]
    if (!(opts.required_arguments({"--in", "--geom", "-o", "-a", "-b"}))) {
        usage(argv[0]);
        return -1;
    }

    auto in_dir = opts.arg("--in");
    auto geoms = opts.args("--geom");
    auto out = opts.arg("-o");
    auto a = opts.arg_as<int>("-a");
    auto b = opts.arg_as<int>("-b");

    auto table = tomo::util::report("Communication volumes", "geometry");
    for (int i = a; i <= b; i *= 2) {
        table.add_column("p = "s + std::to_string(i));
    }

    std::vector<std::thread> threads;
    for (auto i = 0u; i < geoms.size(); ++i) {
        std::cout << in_dir << " " << geoms[i] << "\n";
        threads.emplace_back(compute, in_dir, geoms[i], a, b, std::ref(table),
                             opts.passed("--preview"));
    }

    for (auto& thread : threads) {
        thread.join();
    }

    table.print();

    std::ofstream of(out, std::ios::out);
    table.print(of);

    return 0;
}
