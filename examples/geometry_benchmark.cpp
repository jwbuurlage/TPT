#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <mutex>
#include <thread>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include "fmt/format.h"

#include "tomo.hpp"
#include "util/report.hpp"

std::mutex g_result_mutex;

using T = float;
constexpr tomo::dimension D = 3_D;

void print_node(bulk::binary_tree<bulk::split>::node* node, int depth = 0) {
    for (int i = 0; i < depth; ++i) {
        std::cout << "> ";
    }
    std::cout << "[" << node->value.d << ", " << node->value.a << "]\n";
    if (auto left = node->left.get()) {
        print_node(left, depth + 1);
    }
    if (auto right = node->right.get()) {
        print_node(right, depth + 1);
    }
}

void print_tree(bulk::binary_tree<bulk::split>& tree) {
    print_node(tree.root.get());
}

template <tomo::dimension D, typename T>
void partition_test(std::string name, const tomo::geometry::base<D, T>& g,
                    tomo::util::report& table, tomo::volume<3_D, T> v,
                    int procs, T epsilon) {
    namespace td = tomo::distributed;

    auto part_trivial = td::partition_trivial(g, v, procs);
    auto tree = td::partition_bisection(g, v, procs, epsilon);
    print_tree(tree);
    auto part_bisected = bulk::tree_partitioning<D>(
        tomo::math::vec_to_array<D, int>(v.voxels()), procs, std::move(tree));

    std::cout << "trivial: \n";
    auto overlap_trivial = td::communication_volume<D, T>(g, v, part_trivial);
    std::cout << "bisected: \n";
    auto overlap_bisected = td::communication_volume<D, T>(g, v, part_bisected);

    std::lock_guard<std::mutex> guard(g_result_mutex);
    table.add_row(name);
    table.add_result(name, "trivial", overlap_trivial);
    table.add_result(name, "binary", overlap_bisected);

    T imp = (T)0.0;
    if (overlap_trivial != 0)
        imp = (overlap_trivial - overlap_bisected) / (T)overlap_trivial;

    table.add_result(name, "max_overlap", g.lines() * (procs - 1));
    table.add_result(name, "improvement", fmt::format("{:.1f}%", 100 * imp));
}

int main(int argc, char* argv[]) {
    int k = 32;
    int p = 4;
    T e = (T)0.2;

    po::options_description desc("Allowed arguments");
    desc.add_options()("help,h", "show help message")(
        "size,s", po::value<int>(&k)->default_value(32), "size of the volume")(
        "procs,p", po::value<int>(&p)->default_value(4), "number of procs")(
        "eps,e", po::value<T>(&e)->default_value(e), "load imbalance");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return -1;
    }

    tomo::volume<3_D, T> v(k);

    auto table = tomo::util::report(
        "Detector overlaps for geometric partitioning", "geometry");
    table.add_column("trivial");
    table.add_column("binary");
    table.add_column("max_overlap");
    table.add_column("improvement");

    //    partition_test<D, T>("parallel", tomo::geometry::parallel<3_D, T>(v,
    //    k, k),
    //                         table, v, p, e);
    //
    //    partition_test<D, T>("cone_beam", tomo::geometry::cone_beam<T>(
    //                                          v, k, {1.5, 1.5}, {k, k}, 2.0,
    //                                          2.0),
    //                         table, v, p, e);

    std::vector<std::thread> threads;
    // 'gather results'
    threads.emplace_back(partition_test<D, T>, "parallel",
                         tomo::geometry::parallel<3_D, T>(v, k, k),
                         std::ref(table), v, p, e);

    // threads.emplace_back(
    //     partition_test<D, T>, "dual_parallel",
    //     tomo::geometry::dual_axis_parallel<T>(v, k, {1, 1}, {k, k}),
    //     std::ref(table), v, p, e);

    threads.emplace_back(
        partition_test<D, T>, "dynamic_cone",
        tomo::geometry::dynamic_cone_beam<T>(v, k, {1.5, 1.5}, {k, k}),
        std::ref(table), v, p, e);

    threads.emplace_back(partition_test<D, T>, "helical_cone",
                         tomo::geometry::helical_cone_beam<T>(
                             v, k, {1.5, 1.5}, {k, k}, (T)2.0, (T)2.0, (T)2.0),
                         std::ref(table), v, p, e);

    threads.emplace_back(
        partition_test<D, T>, "cone",
        tomo::geometry::cone_beam<T>(v, k, {1.5, 1.5}, {k, k}, 2.0, 2.0),
        std::ref(table), v, p, e);

    threads.emplace_back(partition_test<D, T>, "laminography",
                         tomo::geometry::laminography<T>(
                             v, k, {1.5, 1.5}, {k, k}, 1.0, 1.0, k / 2, k / 2),
                         std::ref(table), v, p, e);

    threads.emplace_back(
        partition_test<D, T>, "tomo_synthesis",
        tomo::geometry::tomosynthesis<T>(v, k, {1.5, 1.5}, {k, k}),
        std::ref(table), v, p, e);

    for (auto& thread : threads)
        thread.join();

    // 'print result table'
    table.print();

    return 0;
}
