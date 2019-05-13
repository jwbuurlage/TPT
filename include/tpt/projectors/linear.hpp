#pragma once

#include <algorithm>
#include <iostream>

#include "../projector.hpp"
#include "../volume.hpp"

namespace tpt {
namespace dim {

/**
 * This discrete integration method does an interpolation between all the
 * neighbouring voxels of a sample point.
 */
template <dimension D, typename T>
class linear : public base<D, T> {
  public:
    /** Construct the DIM for a given volume. */
    linear(volume<D, T> vol) : base<D, T>(vol) {
        auto max_width = math::max_element<D, T>(vol.voxels());
        this->queue_.reserve((int)(math::sqrt<T>(D) * math::pow(D, 2) * max_width));
    }

    T matrix_value(math::ray<D, T> ray, math::vec<D, int> voxel) {
        (void)ray;
        (void)voxel;
        std::cout << "Warning: `linear::matrix_value` not implemented\n";
        return (T)0;
    }

  private:
    using matrix_iterator = typename base<D, T>::matrix_iterator;

    void reset_(math::line<D, T> line) override {
        auto current_point = line.origin;

        while (math::inside_margin<D, T>(current_point - line.delta,
                                         this->volume_, (T)1.0)) {
            current_point -= line.delta;
        }

        while (math::inside_margin<D, T>(current_point, this->volume_, (T)1.0)) {
            math::interpolate<D, T>(current_point, this->volume_, this->queue_);
            current_point += line.delta;
        }
    }
};

} // namespace dim
} // namespace tpt
