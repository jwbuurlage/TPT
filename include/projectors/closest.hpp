#pragma once

#include "tomo.hpp"

namespace tomo {
namespace dim {

template <typename T>
using closest_iterator =
    typename std::vector<math::matrix_element<T>>::iterator;

/**
 * This discrete integration method attributes the sample point entirely to the
 * closest voxel.
 */
template <dimension D, typename T>
class closest: public base<D, T, closest_iterator<T>> {
  public:
    /** Construct the DIM for a given volume. */
    closest(volume<D> vol)
        : base<D, T, closest_iterator<T>>(vol) {
        auto dims = vol.dimensions();
        auto max_width = *std::max_element(dims.begin(), dims.end());
        queue_.reserve((int)(math::sqrt(D) * max_width));
    }

  private:
    void reset_(line<D, T> line) override {
        queue_.clear();

        auto current_point = line.origin;
        while (math::inside<D, T>(current_point, this->volume_)) {
            std::array<int, D> xs;
            for (int i = 0; i < D; ++i) {
                xs[i] = (int)current_point[i];
            }
            auto index = this->volume_.index(xs);
            if (index >= 0 && index <= this->volume_.cells())
                queue_.push_back({index, (T)1.0});
            current_point += line.delta;
        }

        this->line_ = line;
    }

    closest_iterator<T> begin_() override { return queue_.begin(); }
    closest_iterator<T> end_() override { return queue_.end(); }

    std::vector<math::matrix_element<T>> queue_;
};

} // namespace dim
} // namespace tomo
