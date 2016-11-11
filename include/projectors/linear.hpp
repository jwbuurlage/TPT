#pragma once

#include <algorithm>

#include "../projector.hpp"
#include "../volume.hpp"

namespace tomo {
namespace dim {

/**
 * This discrete integration method does an interpolation between all the
 * neighbouring voxels of a sample point.
 */
template <dimension D, typename T>
class linear
    : public base<D, T,
                  typename std::vector<math::matrix_element<T>>::iterator> {
  public:
    using iterator = typename std::vector<math::matrix_element<T>>::iterator;

    /** Construct the DIM for a given volume. */
    linear(volume<D> vol) : base<D, T, iterator>(vol) {
        auto max_width = math::max_element<D, T>(vol.dimensions());
        queue_.reserve((int)(math::sqrt<T>(D) * math::pow(D, 2) * max_width));
    }

  private:
    void reset_(math::line<D, T> line) override {
        auto current_point = line.origin + (T)0.5 * line.delta;

        while (math::inside<D, T>(current_point, this->volume_)) {
            math::interpolate<D, T>(current_point, this->volume_, queue_);
            current_point += line.delta;
        }

        this->line_ = line;
    }

    void clear_() override { queue_.clear(); }

    iterator begin_() override { return queue_.begin(); }
    iterator end_() override { return queue_.end(); }

    std::vector<math::matrix_element<T>> queue_;
};

} // namespace dim
} // namespace tomo
