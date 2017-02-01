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

template <typename Geometry>
void run(std::string name, Geometry g, tomo::util::report& table,
         tomo::volume<3_D> v, int procs, T epsilon) {
    (void)name;
    (void)g;
    (void)table;
    (void)v;
    (void)procs;
    (void)epsilon;
    // FIXME redo this for the new system
    /*    namespace td = tomo::distributed;
        using Dim = tomo::dim::closest<3_D, T>;

        auto part_trivial = td::partition_trivial<Dim>(g, v, procs);
        auto part_bisected = td::partition_bisection(g, v, procs, epsilon);
        auto overlap_trivial = td::overlap_count<Dim>(g, part_trivial);
        auto overlap_bisected = td::overlap_count<Dim>(g, part_bisected);

        std::lock_guard<std::mutex> guard(g_result_mutex);
        table.add_row(name);
        table.add_result(name, "trivial", overlap_trivial);
        table.add_result(name, "binary", overlap_bisected);

        T imp = (T)0.0;
        if (overlap_trivial != 0)
            imp = (overlap_trivial - overlap_bisected) / (T)overlap_trivial;

        table.add_result(name, "max_overlap", g.lines() * (procs - 1));
        table.add_result(name, "improvement", fmt::format("{:.1f}%", 100 *
       imp)); */
}

int main(int argc, char* argv[]) {
    int k = 64;
    int p = 4;
    T e = (T)0.1;

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

    tomo::volume<3_D> v(k);

    auto table = tomo::util::report(
        "Detector overlaps for geometric partitioning", "geometry");
    table.add_column("trivial");
    table.add_column("binary");
    table.add_column("max_overlap");
    table.add_column("improvement");

    std::vector<std::thread> threads;
    // 'gather results'
    threads.emplace_back(run<tomo::geometry::parallel<3_D, T>>, "parallel",
                         tomo::geometry::parallel<3_D, T>(v, k, k),
                         std::ref(table), v, p, e);

    threads.emplace_back(
        run<tomo::geometry::dual_axis_parallel<T>>, "dual_parallel",
        tomo::geometry::dual_axis_parallel<T>(v, k * 2, (T)1.0, {k, k}),
        std::ref(table), v, p, e);

    threads.emplace_back(
        run<tomo::geometry::dynamic_cone_beam<T>>, "dynamic_cone",
        tomo::geometry::dynamic_cone_beam<T>(v, k, (T)1.0, {k, k}),
        std::ref(table), v, p, e);

    threads.emplace_back(
        run<tomo::geometry::helical_cone_beam<T>>, "helical_cone",
        tomo::geometry::helical_cone_beam<T>(v, k, (T)1.0, {k, k}),
        std::ref(table), v, p, e);

    threads.emplace_back(run<tomo::geometry::cone_beam<T>>, "cone",
                         tomo::geometry::cone_beam<T>(v, k, (T)1.0, {k, k}),
                         std::ref(table), v, p, e);

    threads.emplace_back(run<tomo::geometry::laminography<T>>, "laminography",
                         tomo::geometry::laminography<T>(
                             v, k, (T)1.0, {k, k}, 1.0, 1.0, k / 2, k / 2),
                         std::ref(table), v, p, e);

    threads.emplace_back(run<tomo::geometry::tomosynthesis<T>>,
                         "tomo_synthesis",
                         tomo::geometry::tomosynthesis<T>(v, k, (T)1.0, {k, k}),
                         std::ref(table), v, p, e);

    for (auto& thread : threads)
        thread.join();

    // 'print result table'
    table.print();

    return 0;
}
