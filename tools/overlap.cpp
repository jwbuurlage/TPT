#include <limits>
#include <memory>

#include <filesystem>
namespace fst = std::filesystem;

#include <CLI/CLI.hpp>

#include "tpt/tpt.hpp"

constexpr tpt::dimension D = 3_D;
using T = float;

int main(int argc, char** argv) {
    CLI::App app{"Tool for visualizing the shadow overlap of partitioned "
                 "volumes on the detector"};

    int q = 5;
    int k = (1u << q) * 8u;
    std::vector<std::string> geoms;

    // basic options for partitioning
    app.add_option("-q", q, "log of number of processors to partition for", q);
    app.add_option("-k", k, "number of voxels for stats", k);

    // file input / output
    app.add_option("--geometries", geoms, "the input geometries (.toml files)")
        ->required();

    CLI11_PARSE(app, argc, argv);

    for (auto geom : geoms) {
        auto name = fst::path(geom).stem().string();

        auto problem = tpt::read_configuration<D, T>(geom, k);
        auto v = problem.object_volume;
        auto& g = *problem.acquisition_geometry;

        // plot the overlap for a split in one of the axes
        auto vcorners = v.corners();
        auto corners = std::array<tpt::math::vec3<T>, 8>{};
        std::copy(vcorners.begin(), vcorners.end(), corners.begin());
        auto c = tpt::grcb::split_midway<T>(g, corners, 2);
        auto [vL, vR] = tpt::grcb::split_at(corners, c, 2);
        for (int i = 0; i < q - 2; ++i) {
            auto c = tpt::grcb::split_midway<T>(g, vL, 2);
            auto [vL2, vR2] = tpt::grcb::split_at(vL, c, 2);
            vL = vL2;
            vR = vR2;
        }

        for (auto a : {0, 2}) {
            auto c = tpt::grcb::split_midway<T>(g, vL, a);
            auto [vL2, vR2] = tpt::grcb::split_at(vL, c, a);

            auto A = 0.0f;
            for (auto i = 0; i < g.projection_count(); ++i) {
                A += tpt::grcb::overlap(g.get_projection(i), vL2, vR2);
            }
            std::cout << a << " " << A << "\n";

            for (auto i : {0, 1, 2, 3, 4, 5, 6, 7}) {
                auto pi = g.get_projection(i * (k / 8));

                auto [pL, shadowL] = tpt::grcb::shadow(pi, vL2);
                auto [pR, shadowR] = tpt::grcb::shadow(pi, vR2);

                auto outputname = name + "_overlap_d1_" + std::to_string(i) +
                                  "k3_" + std::to_string(a);

                std::cout << outputname << " "
                          << tpt::grcb::overlap(pi, vL2, vR2) << "\n";

                auto det = tpt::grcb::detector(pi);
                tpt::grcb::plot::overlap<T>(pL, pR, shadowL, shadowR, det,
                                            outputname);
            }
        }

        auto root = tpt::grcb::partition(v, g, q, {});
        tpt::grcb::plot::overlaps(g.get_projection(0), root, corners,
                                  "full_overlap_" + name, 1 << q);

        auto halfvol =
            tpt::volume<3_D, T>({1, 1, 1}, {0.25, 0.25, 0.25}, {0.5, 0.5, 0.5});
        auto hvcorners = halfvol.corners();
        auto hcorners = std::array<tpt::grcb::vec3<T>, 8>{};
        std::copy(hvcorners.begin(), hvcorners.end(), hcorners.begin());
        tpt::grcb::plot::shadow(g.get_projection(2 * k / 3), hcorners,
                                "shadow_" + name);
    }

    return 0;
}
