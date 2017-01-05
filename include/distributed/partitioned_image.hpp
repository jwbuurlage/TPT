#pragma once

#include <algorithm>
#include <chrono>
#include <thread>

#include "bulk/bulk.hpp"

namespace tomo {
namespace distributed {

template <tomo::dimension D>
auto array_to_vec(std::array<int, D> in) {
    tomo::math::vec<D, int> out;
    for (int d = 0; d < D; ++d) {
        out[d] = in[d];
    }
    return out;
}

template <tomo::dimension D, typename T, typename World>
class partitioned_image {
  public:
    partitioned_image(World& world, bulk::partitioning<D, 1>& part)
        : world_(world), part_(part),
          local_volume_(array_to_vec<D>(part_.origin({world.processor_id()})),
                        array_to_vec<D>(part_.local_extent({world.processor_id()}))),
          data_(world_, local_volume_.cells()) {
        clear();
    }

    void clear() { std::fill(data_.begin(), data_.end(), 0); }

    auto& world() const { return world_; }

    tomo::math::vec<D, int> local_size() const {
        return array_to_vec<D>(part_.local_extent({world_.processor_id()}));
    }

    auto global_size() const { return array_to_vec<D>(part_.global_size()); }
    auto local_volume() const { return local_volume_; }

    // get the element with _local index_ idx
    T& operator[](int idx) { return data_[idx]; }

    T& operator()(std::array<int, D> xs) {
        // FIXME this is unnecessarily involved because of the
        // std::array/tomo::math::vec discrepency. Come up with a solution.
        return (*this)[bulk::flatten<D>(
            part_.local_extent({world_.processor_id()}), xs)];
    }

    tomo::image<D, T> gather() {
        auto volume_from_array = [](std::array<int, D> xs) -> tomo::volume<D> {
            tomo::math::vec<D, int> lengths;
            for (int d = 0; d < D; ++d) {
                lengths[d] = xs[d];
            }
            return tomo::volume<D>(lengths);
        };

        auto v = volume_from_array(part_.global_size());
        tomo::image<D, T> result(v);

        int cells = v.cells();
        auto xs = bulk::create_coarray<T>(world_, world_.processor_id() == 0 ? cells : 0);

        int idx = 0;
        for (auto& x : data_) {
            // origin
            // block size
            auto global_pos = bulk::unflatten<D>(
                part_.local_extent({world_.processor_id()}), idx);

            // flatten wrt global size
            for (int d = 0; d < D; ++d) {
                global_pos[d] += part_.origin({world_.processor_id()})[d];
            }

            int global_idx = bulk::flatten<D>(part_.global_size(), global_pos);

            xs(0)[global_idx] = x;

            ++idx;
        }

        world_.sync();

        idx = 0;
        for (auto x : xs) {
            result[idx++] = x;
        }

        return result;
    }

  private:
    World& world_;
    bulk::partitioning<D, 1>& part_;

    tomo::volume<D> local_volume_;
    bulk::coarray<T, World> data_;
};

template <typename Func, typename T, int D, typename World>
void transform_in_place(partitioned_image<D, T, World>& img, Func func) {
    for (auto& elem : img) {
        elem = func(elem);
    }
}

template <typename World, typename Func>
void output_in_turn(World& world, Func func) {
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

template <typename T, typename World>
void partitioned_phantom(partitioned_image<2_D, T, World>& img) {
    auto global_volume =
        tomo::volume<2>(img.global_size()[0], img.global_size()[1]);
    tomo::fill_ellipses_(img, tomo::mshl_ellipses_<T>(), img.local_volume(),
                         global_volume);
}

template <typename T, typename World>
void plot(partitioned_image<2_D, T, World>& img) {
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
