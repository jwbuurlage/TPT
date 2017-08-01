#include "bulk/backends/mpi/mpi.hpp"
#include "bulk/bulk.hpp"

namespace bulk {
using namespace experimental;
}

#include "tomos/tomos.hpp"
using namespace tomo;
using T = float;

struct overlap {
    int projection;
    int target;
    distributed::shadow region;
};

void communicate_overlaps(bulk::world& world,
                          projections<3_D, T>& local_proj_stack,
                          distributed::restricted_geometry<T>& local_geometry,
                          const std::vector<overlap>& overlaps) {
    auto buffer = std::vector<T>();
    auto& proj_data = local_proj_stack.mutable_data();

    auto q = bulk::queue<T[], int, distributed::shadow>(world);
    for (auto ov : overlaps) {
        auto offset = local_proj_stack.offset(ov.projection);
        auto local_region = local_geometry.local_shadow(ov.projection);
        auto w = distributed::shadow{ov.region.min_pt - local_region.min_pt,
                                     ov.region.max_pt - local_region.min_pt};

        buffer.resize(
            math::reduce<2_D>(w.max_pt - w.min_pt + math::vec2<int>{1, 1}));

        auto proj_size = local_geometry.projection_shape(ov.projection);
        int idx = 0;
        for (int j = w.min_pt.y; j <= w.max_pt.y; ++j) {
            for (int i = w.min_pt.x; i <= w.max_pt.x; ++i) {
                buffer[idx++] = proj_data[offset + i + j * proj_size[0]];
            }
        }

        q(ov.target).send_many(buffer, ov.projection, ov.region);
    }

    world.sync();

    for (auto && [ xs, proj, sh ] : q) {
        auto local_region = local_geometry.local_shadow(proj);
        auto w = distributed::shadow{sh.min_pt - local_region.min_pt,
                                     sh.max_pt - local_region.min_pt};

        auto offset = local_proj_stack.offset(proj);
        auto proj_size = local_geometry.projection_shape(proj);
        int idx = 0;
        for (int j = w.min_pt.y; j <= w.max_pt.y; ++j) {
            for (int i = w.min_pt.x; i <= w.max_pt.x; ++i) {
                proj_data[offset + i + j * proj_size[0]] += xs[idx++];
            }
        }
    }
}

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
        auto result = distributed::shadow{};
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
