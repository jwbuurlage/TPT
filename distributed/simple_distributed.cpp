#include <sstream>

#include "bulk/backends/mpi/mpi.hpp"
#include "bulk/bulk.hpp"

namespace bulk {
using namespace experimental;
}

#include "tomos/tomos.hpp"
#include "tomos/util/plotter.hpp"

using namespace tomo;

using T = float;

// Assumes global volume has physical size [0, 1]^3
template <int G>
auto calculate_local_volume(bulk::rectangular_partitioning<3, G>& partitioning,
                            int s) {
    auto to_vec = [](auto array_like) -> math::vec3<T> {
        math::vec3<T> result;
        for (int d = 0; d < 3; ++d) {
            result[d] = array_like[d];
        }
        return result;
    };

    auto voxel_origin = to_vec(partitioning.origin(s));
    auto voxel_size = to_vec(partitioning.local_size(s));
    auto global_voxels = to_vec(partitioning.global_size());

    auto relative_origin = math::vec3<T>(voxel_origin / global_voxels);
    auto relative_size = math::vec3<T>(voxel_size / global_voxels);

    return volume<3_D, T>(voxel_size, relative_origin, relative_size);
}

struct overlap {
    int projection;
    int target;
    distributed::shadow region;
};

template <typename T>
std::vector<overlap>
compute_overlaps(bulk::world& world,
                 const distributed::restricted_geometry<T>& geometry) {

    struct pod_shadow {
        std::array<int, 2> min_pt;
        std::array<int, 2> max_pt;
    };

    auto overlaps = std::vector<overlap>();

    auto s = world.processor_id();
    auto p = world.active_processors();

    // STEP 1: communicate shadows for all projections
    // now have #proj * p shadows
    auto shadows =
        bulk::coarray<pod_shadow>(world, p * geometry.projection_count());

    for (int i = 0; i < geometry.projection_count(); ++i) {
        for (int t = 0; t < p; ++t) {
            auto sh = geometry.local_shadow(i);
            shadows(t)[p * i + s] = {math::vec_to_array<2_D, int>(sh.min_pt),
                                     math::vec_to_array<2_D, int>(sh.max_pt)};
        }
    }
    world.sync();

    auto common_shadow = [](auto s1, auto s2) {
        auto result =
            distributed::shadow{math::array_to_vec<2_D, int>(s1.min_pt),
                                math::array_to_vec<2_D, int>(s1.max_pt)};
        for (int d = 0; d < 2; ++d) {
            result.min_pt[d] = math::max(s1.min_pt[d], s2.min_pt[d]);
            result.max_pt[d] = math::min(s1.max_pt[d], s2.max_pt[d]);
        }
        return result;
    };

    // STEP 2: compute which shadows we have overlap
    //         store these rectangular regions
    for (int i = 0; i < geometry.projection_count(); ++i) {
        for (int t = 0; t < p; ++t) {
            if (t == s) {
                continue;
            }
            auto region = common_shadow(shadows[p * i + s], shadows[p * i + t]);
            if (!distributed::empty(region)) {
                overlaps.push_back({i, t, region});
            }
        }
    }

    // STEP 3: decide who is responsible for summing detector regions
    //         avoid sending overlaps p^2 times
    // TODO

    return overlaps;
}

void run(util::args options) {
    bulk::mpi::environment env;

    // this creates a volume [0, 1]^3, with k voxels in each axis
    auto global_volume = volume<3_D, T>(options.k);
    auto global_geometry =
        geometry::cone_beam<T>(global_volume, options.k, {1.5, 1.5},
                               {options.k, options.k}, 10.0, 2.0);

    auto processors = env.available_processors();
    auto partitioning = bulk::block_partitioning<3, 3>(
        {options.k, options.k, options.k}, {2, 2, 2});

    env.spawn(processors, [&](auto& world) {
        auto local_volume =
            calculate_local_volume(partitioning, world.processor_id());
        (void)world;

        util::ext_plotter<3_D, T> plotter(
            "tcp://localhost:5555", "Distributed restricted geometry test");

        image<3_D, T> phantom(local_volume);
        fill_ellipsoids_(phantom, mshl_ellipsoids_<T>(), local_volume,
                         global_volume);

        auto local_geometry =
            distributed::restricted_geometry(global_geometry, local_volume);
        auto local_proj_stack = projections<3_D, T>(local_geometry);

        // This is now a normal FP
        auto kernel = dim::joseph<3_D, T>(local_volume);
        int line_number = 0;
        for (auto line : local_geometry) {
            for (auto elem : kernel(line)) {
                local_proj_stack[line_number] +=
                    phantom[elem.index] * elem.value;
            }
            ++line_number;
        }

        // ------------------------------------------------------------------
        // overlaps
        auto overlaps = compute_overlaps<T>(world, local_geometry);
        // communicate overlaps where necessary
        auto buffer = std::vector<T>();
        auto& proj_data = local_proj_stack.mutable_data();

        auto q = bulk::queue<T[], int, distributed::shadow>(world);
        for (auto ov : overlaps) {
            auto offset = local_proj_stack.offset(ov.projection);
            auto local_region = local_geometry.local_shadow(ov.projection);
            auto w =
                distributed::shadow{ov.region.min_pt - local_region.min_pt,
                                    ov.region.max_pt - local_region.min_pt};

            buffer.resize(math::reduce<2_D>(w.max_pt - w.min_pt));

            auto proj_size = local_geometry.projection_shape(ov.projection);
            int idx = 0;
            for (int i = w.min_pt.x; i < w.max_pt.x; ++i) {
                for (int j = w.min_pt.y; j < w.max_pt.y; ++j) {
                    buffer[idx++] = proj_data[offset + i + j * proj_size[0]];
                }
            }

            q(ov.target).send_many(buffer, ov.projection, ov.region);
        }

        world.sync();

        for (auto&& [xs, proj, sh] : q) {
            auto local_region = local_geometry.local_shadow(proj);
            auto w = distributed::shadow{sh.min_pt - local_region.min_pt,
                                         sh.max_pt - local_region.min_pt};
            auto offset = local_proj_stack.offset(proj);
            auto proj_size = local_geometry.projection_shape(proj);
            int idx = 0;
            for (int i = w.min_pt.x; i < w.max_pt.x; ++i) {
                for (int j = w.min_pt.y; j < w.max_pt.y; ++j) {
                    proj_data[offset + i + j * proj_size[0]] = xs[idx++];
                }
            }
        }
        // ------------------------------------------------------------------

        plotter.plot(phantom);
        plotter.send_projection_data(local_geometry, local_proj_stack,
                                     local_volume);
    });
}

int main(int argc, char* argv[]) {
    auto options = util::args(argc, argv);
    run(options);

    return 0;
}
