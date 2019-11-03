#pragma once

#include <limits>

#include "../geometry.hpp"
#include "../math.hpp"
#include "../math/stringify.hpp"
#include "../volume.hpp"
#include "trajectory.hpp"

namespace tpt {
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
                      math::vec<2_D, int> detector_shape,
                      math::vec<3_D, T> source_position,
                      math::vec<3_D, T> detector_position, T rotations,
                      std::array<math::vec<3_D, T>, 2> detector_tilt,
                      std::vector<T> angles = {})
        : trajectory<3_D, T>(volume, projection_count, detector_size,
                             detector_shape),
          source_position_(source_position),
          detector_position_(detector_position), rotations_(rotations),
          detector_tilt_(detector_tilt), angles_(std::move(angles)) {
        assert(projection_count > 1);
        if (angles_.empty()) {
            T delta =
                (rotations_ * (T)2.0 * math::pi<T>) / (projection_count - 1);
            for (int i = 0; i < projection_count; ++i) {
                angles_.push_back(delta * i);
            }
        }
        assert((int)angles_.size() == projection_count);
    }

    /** Construct the geometry with a given number of lines. */
    helical_cone_beam(volume<3_D, T> volume, int projection_count,
                      math::vec<2_D, T> detector_size,
                      math::vec<2_D, int> detector_shape, T rotations,
                      T source_to_center, T detector_to_center)
        : helical_cone_beam(
              volume, projection_count, detector_size, detector_shape,
              math::volume_center(volume) -
                  source_to_center * math::standard_basis<3_D, T>(0),
              math::volume_center(volume) +
                  detector_to_center * math::standard_basis<3_D, T>(0),
              rotations,
              {math::standard_basis<3_D, T>(1),
               math::standard_basis<3_D, T>(2)}) {}

    math::vec<3_D, T> source_location(int projection) const override final {
        return transform_location_(source_position_, projection) +
               height_(projection);
    }

    math::vec<3_D, T> detector_location(int projection) const override final {
        return transform_location_(detector_position_, projection) +
               height_(projection);
    }

    std::array<math::vec<3_D, T>, 2>
    detector_tilt(int projection) const override final {
        return {apply_rotation_(detector_tilt_[0], projection),
                apply_rotation_(detector_tilt_[1], projection)};
    }

  private:
    math::vec<3_D, T> source_position_;
    math::vec<3_D, T> detector_position_;

    T rotations_;

    std::array<math::vec<3_D, T>, 2> detector_tilt_;
    std::vector<T> angles_;

    inline math::vec<3_D, T> transform_location_(math::vec<3_D, T> location,
                                                 int projection) const {
        return apply_rotation_(location - math::volume_center(this->volume_),
                               projection) +
               math::volume_center(this->volume_);
    }

    inline math::vec<3_D, T> apply_rotation_(math::vec<3_D, T> location,
                                             int projection) const {
        return math::rotate(location, math::standard_basis<3_D, T>(2),
                            angles_[projection]);
    }

    inline math::vec<3_D, T> height_(int projection) const {
        return ((((T)projection / (this->projection_count_ - 1)) - (T)0.5) *
                this->volume_[2]) *
               math::standard_basis<3_D, T>(2);
    }
};

} // namespace geometry
} // namespace tpt
