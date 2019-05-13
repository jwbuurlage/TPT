#pragma once

#include "../common.hpp"
#include "../math.hpp"
#include "../projector.hpp"

namespace tpt {
namespace dim {

/**
 * This discrete integration method attributes the sample point entirely to the
 * closest voxel.
 */
template <dimension D, typename T>
class closest : public base<D, T> {
  public:
    /** Construct the DIM for a given volume. */
    closest(volume<D, T> vol) : base<D, T>(vol) {
        auto max_width = math::max_element<D, T>(vol.voxels());
        this->queue_.reserve((int)(math::sqrt(D) * max_width));
    }


    T matrix_value(math::ray<D, T> ray, math::vec<D, int> voxel) {
        (void)ray;
        (void)voxel;
        std::cout << "Warning: `closest::matrix_value` not implemented\n";
        return (T)0;
    }

  private:
    using matrix_iterator = typename base<D, T>::matrix_iterator;

    void reset_(math::line<D, T> line) override {
        auto current_point = line.origin + (T)0.5 * line.delta;
        while (math::inside<D, T>(current_point, this->volume_)) {
            // convert to vector of integers
            auto index = (uint64_t)this->volume_.index(
                math::vec<D, int>(current_point));
            if (index >= 0 && index < this->volume_.cells()) {
                this->queue_.push_back({(int)index, (T)1.0});
            }
            current_point += line.delta;
        }
    }
};

} // namespace dim
} // namespace tpt
