#pragma once

#include "bulk/bulk.hpp"
namespace bulk {
using namespace experimental;
}

#include "../math.hpp"
#include "../projectors/closest.hpp"

namespace tomo {
namespace distributed {

template <dimension D, typename T>
long long communication_volume(const geometry::base<D, T>& geometry,
                               volume<D, T> object_volume,
                               bulk::partitioning<D>& partitioning) {
    long long result = 0;

    auto integrator = tomo::dim::closest<D, T>(object_volume);
    auto voxels = math::vec_to_array<D, int>(object_volume.voxels());

    std::set<int> owners;
    for (auto [idx, line] : geometry) {
        (void)idx;
        owners.clear();
        for (auto elem : integrator(line)) {
            auto voxel_idx = bulk::util::unflatten<D>(voxels, elem.index);
            auto owner = partitioning.owner(voxel_idx);
            owners.insert(owner);
        }
        result += math::max(0, (int)owners.size() - 1);
    }

    return result;
}

template <dimension D, typename T, dimension G>
T load_imbalance(tomo::volume<D, T> v,
                 bulk::rectangular_partitioning<D, G>& partitioning,
                 const geometry::base<D, T>& geometry) {
    auto p =
        math::product<G, int>(math::array_to_vec<G, int>(partitioning.grid()));
    auto ws = voxel_weights(geometry, v);
    auto avg_weight =
        weight({0, 0, 0},
               math::array_to_vec<D, int>(partitioning.global_size()) -
                   math::vec<D, int>{1},
               ws) /
        p;

    // for each owner, get bounds
    auto eps = (T)0;
    for (int s = 0; s < p; ++s) {
        auto base = math::array_to_vec<D, int>(partitioning.origin(s));
        auto end = base +
                   math::array_to_vec<D, int>(partitioning.local_size(s)) -
                   math::vec<D, int>{1};
        auto w = weight(base, end, ws);
        eps = math::max(w / avg_weight - (T)1, eps);
    }

    return eps;
}

template <dimension D, typename T, dimension G>
long long regularizer_volume(volume<D, T> object_volume,
                   bulk::rectangular_partitioning<D, G>& partitioning) {
    long long result = 0;
    auto p =
        math::product<G, int>(math::array_to_vec<G, int>(partitioning.grid()));

    auto boundary_voxels = [](auto& xs) {
        return 2 * (xs[0] * xs[1] + (xs[0] - 1) * xs[2] + (xs[1] - 1) * (xs[2] - 1))
    };

    for (int s = 0; s < p; ++s) {
        auto xs = math::array_to_vec<D, int>(partitioning.local_size(s));
        result += boundary_voxels(xs);
    }

    auto ys = math::vec_to_array<D, int>(object_volume.voxels());
    result -= boundary_voxels(ys);

    return result;
}

} // namespace distributed
} // namespace tomo
