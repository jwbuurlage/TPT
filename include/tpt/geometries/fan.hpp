#pragma once

#include <limits>

#include "../geometry.hpp"
#include "../math.hpp"
#include "../volume.hpp"
#include "trajectory.hpp"

namespace tpt {
namespace geometry {

/**
 * A circular fan-beam geometry
 *
 * \tparam T the scalar type to use
 */
template <typename T>
class fan_beam : public trajectory<2_D, T> {
  public:
    /** Construct the geometry with a given number of lines. */
    fan_beam(volume<2_D, T> volume, int projection_count,
             math::vec<1_D, T> detector_size,
             math::vec<1_D, int> detector_shape,
             math::vec<2_D, T> source_position,
             math::vec<2_D, T> detector_position,
             std::array<math::vec<2_D, T>, 1> detector_tilt,
             std::vector<T> angles = {})
        : trajectory<2_D, T>(volume, projection_count, detector_size,
                             detector_shape),
          source_position_(source_position),
          detector_position_(detector_position), detector_tilt_(detector_tilt),
          angles_(std::move(angles)) {
        if (angles_.empty()) {
            T delta = ((T)2.0 * math::pi<T>) / projection_count;
            for (int i = 0; i < projection_count; ++i) {
                angles_.push_back(delta * i);
            }
        }
        assert((int)angles_.size() == projection_count);
    }

    /** Construct a fan_beam geometry using many defaults. */
    fan_beam(volume<2_D, T> volume, int projection_count,
             math::vec<1_D, T> detector_size,
             math::vec<1_D, int> detector_shape, T source_to_center,
             T detector_to_center)
        : fan_beam(volume, projection_count, detector_size, detector_shape,
                   math::volume_center(volume) -
                       source_to_center * math::standard_basis<2_D, T>(0),
                   math::volume_center(volume) +
                       detector_to_center * math::standard_basis<2_D, T>(0),
                   {math::standard_basis<2_D, T>(1)}) {}

    math::vec<2_D, T> source_location(int projection) const override final {
        return transform_location_(source_position_, projection);
    }

    math::vec<2_D, T> detector_location(int projection) const override final {
        return transform_location_(detector_position_, projection);
    }

    std::array<math::vec<2_D, T>, 1>
    detector_tilt(int projection) const override final {
        return {apply_rotation_(detector_tilt_[0], projection)};
    }

  private:
    math::vec<2_D, T> source_position_;
    math::vec<2_D, T> detector_position_;
    std::array<math::vec<2_D, T>, 1> detector_tilt_;
    std::vector<T> angles_;

    inline math::vec<2_D, T> transform_location_(math::vec<2_D, T> location,
                                                 int projection) const {
        return apply_rotation_(location - math::volume_center(this->volume_),
                               projection) +
               math::volume_center(this->volume_);
    }

    inline math::vec<2_D, T> apply_rotation_(math::vec<2_D, T> location,
                                             int projection) const {
        return math::rotate(location, math::standard_basis<2_D, T>(1),
                            angles_[projection]);
    }
};

} // namespace geometry
} // namespace tpt
