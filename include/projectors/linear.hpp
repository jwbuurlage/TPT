#pragma once

#include <algorithm>

#include "../projector.hpp"
#include "../volume.hpp"

namespace tomo {

/**
 * This discrete integration method does an interpolation between all the
 * neighbouring voxels of a sample point.
 */
template <dimension D, typename T>
class linear_projector
    : public projector<
          D, T, typename std::vector<math::matrix_element<T>>::iterator> {
  public:
    using iterator = typename std::vector<math::matrix_element<T>>::iterator;

    /** Construct the DIM for a given volume. */
    linear_projector(volume<D> vol) : projector<D, T, iterator>(vol) {
        auto dims = this->volume_.dimensions();
        auto max_width = *std::max_element(dims.begin(), dims.end());
        queue_.reserve((int)(math::sqrt<T>(D) * math::pow(D, 2) * max_width));
    }

  private:
    void reset_(line<D, T> line) override {
        queue_.clear();
        auto current_point = line.origin;

        while (math::inside<D, T>(current_point, this->volume_)) {
            math::interpolate<D, T>(current_point, this->volume_, queue_);
            current_point += line.delta;
        }

        this->line_ = line;
    }

    iterator begin_() override { return queue_.begin(); }
    iterator end_() override { return queue_.end(); }

    std::vector<math::matrix_element<T>> queue_;
};

} // namespace tomo
