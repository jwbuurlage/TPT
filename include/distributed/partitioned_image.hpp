#pragma once

#include <algorithm>
#include <chrono>
#include <thread>

#include "../phantoms.hpp"
#include "bulk/bulk.hpp"

namespace tomo {
namespace distributed {

/** Partitioned image partitions the _voxels_ */
template <tomo::dimension D, tomo::dimension G, typename T>
class partitioned_image {
  public:
    partitioned_image(bulk::world& world,
                      bulk::rectangular_partitioning<D, G>& part,
                      tomo::volume<D, T> global_vol)
        : world_(world), part_(part),
          data_(world_, part.local_count(world.processor_id())),
          local_volume_(0) {
        auto voxels =
            math::array_to_vec<D, int>(part_.local_size(world.processor_id()));
        auto voxel_origin =
            math::array_to_vec<D>(part_.origin({world.processor_id()}));

        auto relative_origin = math::vec<D, T>(voxel_origin) /
                               math::vec<D, T>(global_vol.voxels());
        auto origin = global_vol.origin() +
                      relative_origin * global_vol.physical_lengths();

        auto lengths =
            (math::vec<D, T>(voxels) / math::vec<D, T>(global_vol.voxels())) *
            global_vol.physical_lengths();

        local_volume_ = tomo::volume<D, T>(voxels, origin, lengths);
        clear();
    }

    void clear() { std::fill(data_.begin(), data_.end(), 0); }

    auto& world() const { return world_; }

    tomo::math::vec<D, int> local_size() const {
        return math::array_to_vec<D>(part_.local_size(world_.processor_id()));
    }

    auto global_size() const {
        return math::array_to_vec<D>(part_.global_size());
    }
    auto local_volume() const { return local_volume_; }

    // get the element with _local index_ idx
    T& operator[](int idx) { return data_[idx]; }

    T& operator()(std::array<int, D> xs) {
        // FIXME this is unnecessarily involved because of the
        // std::array/tomo::math::vec discrepency. Come up with a solution.
        return (*this)[bulk::flatten<D>(part_.local_size(world_.processor_id()),
                                        xs)];
    }

    tomo::image<D, T> gather() {
        auto volume_from_array =
            [](std::array<int, D> xs) -> tomo::volume<D, T> {
            tomo::math::vec<D, int> lengths;
            for (int d = 0; d < D; ++d) {
                lengths[d] = xs[d];
            }
            return tomo::volume<D, T>(lengths);
        };

        auto v = volume_from_array(part_.global_size());
        tomo::image<D, T> result(v);

        int cells = v.cells();
        auto xs =
            bulk::coarray<T>(world_, world_.processor_id() == 0 ? cells : 0);


        int local_idx = 0;
        for (auto& x : data_) {
            // origin
            // block size
            auto global_pos = bulk::unflatten<D>(
                part_.local_size(world_.processor_id()), local_idx);

            // flatten wrt global size
            for (int d = 0; d < D; ++d) {
                global_pos[d] += part_.origin({world_.processor_id()})[d];
            }

            int global_idx = bulk::flatten<D>(part_.global_size(), global_pos);

            xs(0)[global_idx] = x;

            ++local_idx;
        }

        world_.sync();

        int idx = 0;
        for (auto x : xs) {
            result[idx++] = x;
        }

        return result;
    }

  private:
    bulk::world& world_;
    bulk::rectangular_partitioning<D, G>& part_;
    bulk::coarray<T> data_;

    tomo::volume<D, T> local_volume_;
};

template <typename Func, typename T, int D, int G>
void transform_in_place(partitioned_image<D, G, T>& img, Func func) {
    for (auto& elem : img) {
        elem = func(elem);
    }
}

template <typename Func>
void output_in_turn(bulk::world& world, Func func) {
    // easiest is to copy the local elements into a tomo::image, then output the
    // pid, origin, extent, and finally plot it.
    // we ensure correct oder here
    for (int t = 0; t < world.active_processors(); ++t) {
        if (world.processor_id() == t) {
            func();
        }

        std::chrono::milliseconds timespan(100);
        std::this_thread::sleep_for(timespan);

        world.sync();
    }
}

template <typename T, int G>
void partitioned_phantom(partitioned_image<2_D, G, T>& img) {
    auto global_volume = tomo::volume<2_D, T>(img.global_size());
    tomo::fill_ellipses_(img, tomo::mshl_ellipses_<T>(), img.local_volume(),
                         global_volume);
}

template <typename T, int G>
void partitioned_phantom(partitioned_image<3_D, G, T>& img) {
    auto global_volume = tomo::volume<3_D, T>(img.global_size());
    tomo::fill_ellipsoids_(img, tomo::mshl_ellipsoids_<T>(), img.local_volume(),
                           global_volume);
}

template <typename T, int G>
void plot(partitioned_image<2_D, G, T>& img) {
    auto& world = img.world();
    auto volume = img.local_volume();
    auto cells = volume.cells();

    T max = std::numeric_limits<T>::min();

    auto local_image = tomo::image<2_D, T>(img.local_volume());
    for (int i = 0; i < cells; ++i) {
        local_image[i] = img[i];
        if (max < img[i]) {
            max = img[i];
        }
    }

    // exchange global maximum
    auto maxs = bulk::gather_all(world, max);
    max = *std::max_element(maxs.begin(), maxs.end());

    output_in_turn(world, [&]() {
        world.log("(%i x %i) + (%i, %i) \n", img.local_size()[0],
                  img.local_size()[1], volume.origin()[0], volume.origin()[1]);
        tomo::ascii_plot<T>(local_image, max);
    });
}

} // namespace distributed
} // namespace tomo
