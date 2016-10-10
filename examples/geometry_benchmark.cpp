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

template <typename Geometry>
void run(std::string name, Geometry g, tomo::util::report& table,
         tomo::volume<3_D> v) {
    table.add_row(name);

    auto part_trivial = tomo::distributed::partition_trivial(g, v, 4);
    table.add_result(name, "trivial", overlap_count(g, part_trivial));

    auto part_bisected = tomo::distributed::partition_bisection(g, v, 4);
    table.add_result(name, "binary", overlap_count(g, part_bisected));
}

int main(int argc, char* argv[]) {
    using T = float;

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
    // 5) [ ] make the recursive bijection partitioner and compare
    // 6) [ ] add cli flags

    int k = 16;
    tomo::volume<3_D> v(k);

    auto table = tomo::util::report(
        "Detector overlaps for geometric partitioning", "geometry");
    table.add_column("trivial");
    table.add_column("binary");

    // 'gather results'
    run("parallel", tomo::geometry::parallel<3_D, T>(k, k, v), table, v);
    run("dual_parallel",
        tomo::geometry::dual_axis_parallel<T>(v, k, (T)1.0, {k, k}), table, v);
    run("cone", tomo::geometry::cone_beam<T>(v, k, (T)1.0, {k, k}), table, v);
    run("laminography", tomo::geometry::laminography<T>(v, k, (T)1.0, {k, k},
                                                       1.0, 1.0, k / 2, k / 2),
        table, v);

    // 'print result table'
    table.print();

    return 0;
}
