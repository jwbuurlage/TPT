#pragma once

#include <limits>

#include "../geometry.hpp"
#include "../math.hpp"
#include "../volume.hpp"
#include "trajectory.hpp"

namespace tomo {
namespace geometry {

/**
 * A parallel geometry over two axes
 *
 * \tparam D the dimension of the volume.
 * \tparam T the scalar type to use
 */
template <typename T>
class dual_axis_parallel : public trajectory<3_D, T> {
  public:
    dual_axis_parallel(
        volume<3_D, T> volume, int steps, T detector_spacing = (T)1,
        math::vec<2_D, int> detector_size = math::vec<2_D, int>{1})
        : trajectory<3_D, T>(volume,
                             2 * steps * math::reduce<2_D>(detector_size),
                             detector_spacing, math::vec<2_D, int>{1}),
          actual_detector_size_(detector_size) {
        detector_count_copy_ = math::reduce<2_D>(detector_size);
        actual_steps_ = steps;

        assert(detector_count_copy_ > 0);
        assert(actual_steps_ > 0);
    }

    math::vec<3_D, T> source_location(int step) const override final {
        return transform_location_(
                   math::volume_center(this->volume_) -
                       ((T) this->volume_[0] * math::standard_basis<3_D, T>(0)),
                   step_actual_(step)) +
               this->detector_offset_(step, step % detector_count_copy_,
                                      actual_detector_size_);
    }

    math::vec<3_D, T> detector_location(int step) const override final {
        return transform_location_(
                   math::volume_center(this->volume_) +
                       ((T) this->volume_[0] * math::standard_basis<3_D, T>(0)),
                   step_actual_(step)) +
               this->detector_offset_(step, step % detector_count_copy_,
                                      actual_detector_size_);
    }

    std::array<math::vec<3_D, T>, 2>
    detector_tilt(int step) const override final {
        return {apply_rotation_(math::standard_basis<3_D, T>(1),
                                step_actual_(step)),
                apply_rotation_(math::standard_basis<3_D, T>(2),
                                step_actual_(step))};
    }

  private:
    inline math::vec<3_D, T> transform_location_(math::vec<3_D, T> location,
                                                 int step) const {
        return apply_rotation_(location - math::volume_center(this->volume_), step) +
                math::volume_center(this->volume_);
    }

    inline math::vec<3_D, T> apply_rotation_(math::vec<3_D, T> location,
                                             int step) const {
        static auto axes = std::array<math::vec<3_D, T>, 2>{
            math::standard_basis<3_D, T>(1), math::standard_basis<3_D, T>(2)};
        auto axis = axes[step / actual_steps_];

        T angle_step = math::pi<T> / actual_steps_;

        return math::rotate(location, axis, angle_step * (step % actual_steps_));
    }

    int step_actual_(int step) const { return step / detector_count_copy_; }

    int detector_count_copy_;
    int actual_steps_;
    math::vec<2_D, int> actual_detector_size_;
};

} // namespace geometry
} // namespace tomo
