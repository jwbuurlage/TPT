#pragma once

#include <limits>

#include "../geometry.hpp"
#include "../math.hpp"
#include "../volume.hpp"
#include "trajectory.hpp"

namespace tomo {
namespace geometry {

/**
 * A helical cone-beam geometry
 *
 * \tparam D the dimension of the volume.
 * \tparam T the scalar type to use
 */
template <typename T>
class helical_cone_beam : public trajectory<3_D, T> {
  public:
    /** Construct the geometry with a given number of lines. */
    helical_cone_beam(
        volume<3_D, T> volume, int steps, T detector_spacing = (T)1.0,
        math::vec<2_D, int> detector_size = math::vec<2_D, int>{1},
        T relative_source_distance = (T)1.0,
        T relative_detector_distance = (T)1.0)
        : trajectory<3_D, T>(volume, steps, detector_spacing, detector_size),
          relative_source_distance_(relative_source_distance),
          relative_detector_distance_(relative_detector_distance) {}

    math::vec<3_D, T> source_location(int step) const override final {
        return transform_location_(math::volume_center(this->volume_) -
                                       (relative_source_distance_ *
                                        this->volume_[0] *
                                        math::standard_basis<3_D, T>(0)),
                                   step) +
               height_(step);
    }

    math::vec<3_D, T> detector_location(int step) const override final {
        return transform_location_(math::volume_center(this->volume_) +
                                       (relative_detector_distance_ *
                                        this->volume_[0] *
                                        math::standard_basis<3_D, T>(0)),
                                   step) +
               height_(step);
    }

    std::array<math::vec<3_D, T>, 2>
    detector_tilt(int step) const override final {
        return {apply_rotation_(math::standard_basis<3_D, T>(1), step),
                apply_rotation_(math::standard_basis<3_D, T>(2), step)};
    }

    T& relative_source_distance() { return relative_source_distance_; }
    T& relative_detector_distance() { return relative_detector_distance_; }

  private:
    T relative_source_distance_;
    T relative_detector_distance_;

    inline math::vec<3_D, T> transform_location_(math::vec<3_D, T> location,
                                                 int step) const {
        return apply_rotation_(location - math::volume_center(this->volume_), step) +
               math::volume_center(this->volume_);
    }

    inline math::vec<3_D, T> apply_rotation_(math::vec<3_D, T> location,
                                             int step) const {
        static auto axis = math::standard_basis<3_D, T>(1);

        T angle_step = (T)2.0 * math::pi<T> / this->steps_;

        return math::rotate(location, axis, angle_step * step);
    }

    inline math::vec<3_D, T> height_(int step) const {
        return ((((T)step / this->steps_) - (T)0.5) * this->volume_[2]) *
               math::standard_basis<3_D, T>(1);
    }
};

} // namespace geometry
} // namespace tomo
