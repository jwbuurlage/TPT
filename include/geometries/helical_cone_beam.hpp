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
    helical_cone_beam(volume<3_D, T> volume, int projection_count,
                      math::vec<2_D, T> detector_size,
                      math::vec<2_D, int> detector_shape, T rotations,
                      T relative_source_distance = 1.0,
                      T relative_detector_distance = 1.0)
        : trajectory<3_D, T>(volume, projection_count, detector_size,
                             detector_shape),
          rotations_(rotations),
          relative_source_distance_(relative_source_distance),
          relative_detector_distance_(relative_detector_distance) {}

    math::vec<3_D, T> source_location(int projection) const override final {
        return transform_location_(math::volume_center(this->volume_) -
                                       (relative_source_distance_ *
                                        this->volume_[0] *
                                        math::standard_basis<3_D, T>(0)),
                                   projection) +
               height_(projection);
    }

    math::vec<3_D, T> detector_location(int projection) const override final {
        return transform_location_(math::volume_center(this->volume_) +
                                       (relative_detector_distance_ *
                                        this->volume_[0] *
                                        math::standard_basis<3_D, T>(0)),
                                   projection) +
               height_(projection);
    }

    std::array<math::vec<3_D, T>, 2>
    detector_tilt(int projection) const override final {
        return {apply_rotation_(math::standard_basis<3_D, T>(1), projection),
                apply_rotation_(math::standard_basis<3_D, T>(2), projection)};
    }

    T& relative_source_distance() { return relative_source_distance_; }
    T& relative_detector_distance() { return relative_detector_distance_; }

  private:
    T rotations_;
    T relative_source_distance_;
    T relative_detector_distance_;

    inline math::vec<3_D, T> transform_location_(math::vec<3_D, T> location,
                                                 int projection) const {
        return apply_rotation_(location - math::volume_center(this->volume_),
                               projection) +
               math::volume_center(this->volume_);
    }

    inline math::vec<3_D, T> apply_rotation_(math::vec<3_D, T> location,
                                             int projection) const {
        static auto axis = math::standard_basis<3_D, T>(1);

        T angle_projection = rotations_ * (T)2.0 * math::pi<T> / this->projection_count_;

        return math::rotate(location, axis, angle_projection * projection);
    }

    inline math::vec<3_D, T> height_(int projection) const {
        return ((((T)projection / this->projection_count_) - (T)0.5) *
                this->volume_[2]) *
               math::standard_basis<3_D, T>(1);
    }
};

} // namespace geometry
} // namespace tomo
