#pragma once

#include <limits>

#include "partitioning_statistics.hpp"

namespace tpt {
namespace distributed {

template <dimension D, typename T>
bulk::block_partitioning<D, 1> partition_trivial(const geometry::base<D, T>& g,
                                                 volume<D, T> v,
                                                 int processors) {
    bulk::block_partitioning<D, 1> result(
        math::vec_to_array<D, int>(v.voxels()), {processors}, {0});

    int best_communication_volume = std::numeric_limits<int>::max();
    auto voxels = math::vec_to_array<D, int>(v.voxels());
    for (int axis = 0; axis < D; ++axis) {
        auto part =
            bulk::block_partitioning<D, 1>(voxels, {processors}, {axis});
        auto com_vol = communication_volume(g, v, part);
        if (com_vol < best_communication_volume) {
            best_communication_volume = com_vol;
            result = part;
        }
    }

    return result;
}

} // namespace distributed
} // namespace tpt
