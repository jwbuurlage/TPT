#pragma once

#include "../common.hpp"
#include "../math.hpp"
#include "../projector.hpp"

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
class closest : public base<D, T, closest_iterator<T>> {
  public:
    /** Construct the DIM for a given volume. */
    closest(volume<D> vol) : base<D, T, closest_iterator<T>>(vol) {
        auto max_width = math::max_element<D, T>(vol.dimensions());
        queue_.reserve((int)(math::sqrt(D) * max_width));
    }

  private:
    void reset_(math::line<D, T> line) override {
        auto current_point = line.origin;
        while (math::inside<D, T>(current_point, this->volume_)) {
            // implicitely convert to vector of integers
            auto index = this->volume_.index(math::vec<D, int>(current_point) -
                                             this->volume_.origin());
            if (index >= 0 && index <= this->volume_.cells()) {
                queue_.push_back({index, (T)1.0});
            }
            current_point += line.delta;
        }

        this->line_ = line;
    }

    void clear_() override { queue_.clear(); }

    closest_iterator<T> begin_() override { return queue_.begin(); }
    closest_iterator<T> end_() override { return queue_.end(); }

    std::vector<math::matrix_element<T>> queue_;
};

} // namespace dim
} // namespace tomo
