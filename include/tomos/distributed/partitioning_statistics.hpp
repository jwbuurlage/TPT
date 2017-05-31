#pragma once

#include "bulk/bulk.hpp"
namespace bulk {
using namespace experimental;
}

#include "../projectors/closest.hpp"

namespace tomo {
namespace distributed {

template <dimension D, typename T>
int communication_volume(const geometry::base<D, T>& geometry,
                         volume<D, T> object_volume,
                         bulk::partitioning<D>& partitioning) {
    int result = 0;

    // TODO comvol is now wrt closest, how do we want to do this
    // 'honest' thing is maybe wrt line kernel
    auto integrator = tomo::dim::closest<D, T>(object_volume);
    auto voxels = math::vec_to_array<D, int>(object_volume.voxels());

    for (auto line : geometry) {
        int count = 0;
        int prev = -1;
        for (auto elem : integrator(line)) {
            auto owner = partitioning.owner(
                bulk::util::unflatten<D>(voxels, elem.index));
            if (owner != prev) {
                count++;
                prev = owner;
            }
        }
        result += math::max(0, count - 1);
    }

    return result;
}

template <dimension D, typename T>
T load_imbalance(tomo::volume<D, T> v, bulk::partitioning<D>& partitioning) {
    (void)v;
    (void)partitioning;
    return (T)0.0;
}

} // namespace distributed
} // namespace tomo
