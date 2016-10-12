#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include "fmt/format.h"

#include <boost/hana.hpp>
namespace hana = boost::hana;

#include "tomo.hpp"
#include "util/report.hpp"

using T = float;

template <typename Geometry>
void run(std::string name, Geometry g, tomo::util::report& table,
         tomo::volume<3_D> v, int procs) {
    namespace td = tomo::distributed;
    using Dim = tomo::dim::closest<3_D, T>;

    table.add_row(name);

    auto part_trivial = td::partition_trivial<Dim>(g, v, procs);
    table.add_result(name, "trivial", td::overlap_count<Dim>(g, part_trivial));

    auto part_bisected = td::partition_bisection(g, v, procs);
    table.add_result(name, "binary", td::overlap_count<Dim>(g, part_bisected));

    table.add_result(name, "max_overlap", g.lines() * (procs - 1));
}

int main(int argc, char* argv[]) {
    (void)argc;
    (void)argv;
    // we need simple mechanism for computing 'size of overlap', and see what
    // trivial
    // simple partitioners yield. Also, visualize result, think of 'distributed
    // volume' object.
    // 1) [x] compute the detector extent, or which lines intersect the volume
    // 2) [x] sum over all overlaps
    // 3) [x] support a distribution of the volume
    // 4) [x] make a trivial partition[er/ing]
    // 5) [x] make the recursive bijection partitioner and compare
    // 6) [x] add cli flags

    int k = 64;
    int p = 4;

    po::options_description desc("Allowed arguments");
    desc.add_options()("help,h", "show help message")(
        "size,s", po::value<int>(&k)->default_value(32), "size of the volume")(
        "procs,p", po::value<int>(&p)->default_value(4), "number of procs");

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

    // 'gather results'
    run("parallel", tomo::geometry::parallel<3_D, T>(k, k, v), table, v, p);
    run("dual_parallel",
        tomo::geometry::dual_axis_parallel<T>(v, k * 2, (T)1.0, {k, k}), table,
        v, p);
    run("dynamic_cone", tomo::geometry::cone_beam<T>(v, k, (T)1.0, {k, k}),
        table, v, p);
    run("helical_cone",
        tomo::geometry::helical_cone_beam<T>(v, k, (T)1.0, {k, k}), table, v, p);
    run("cone", tomo::geometry::cone_beam<T>(v, k, (T)1.0, {k, k}), table, v, p);
    run("laminography", tomo::geometry::laminography<T>(v, k, (T)1.0, {k, k},
                                                        1.0, 1.0, k / 2, k / 2),
        table, v, p);
    run("tomo_synthesis",
        tomo::geometry::tomosynthesis<T>(v, k, (T)1.0, {k, k}), table, v, p);

    // 'print result table'
    table.print();

    return 0;
}
