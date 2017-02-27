#pragma once

#include <limits>

#include "../geometry.hpp"
#include "../math.hpp"
#include "../volume.hpp"
#include "trajectory.hpp"

namespace tomo {
namespace geometry {

/**
 * A circular cone-beam geometry
 *
 * \tparam D the dimension of the volume.
 * \tparam T the scalar type to use
 */
template <typename T>
class dynamic_cone_beam : public trajectory<3_D, T> {
  public:
    /** Construct the geometry with a given number of lines. */
    dynamic_cone_beam(volume<3_D, T> volume, int projection_count,
              math::vec<2_D, T> detector_size,
              math::vec<2_D, int> detector_shape)
        : trajectory<3_D, T>(volume, projection_count, detector_size,
                             detector_shape) {}

    math::vec<3_D, T> source_location(int projection) const override final {
        return transform_location_(image_center_() -
                                       (relative_source_distance_(projection) *
                                        this->volume_[0] *
                                        math::standard_basis<3_D, T>(0)),
                                   projection);
    }

    math::vec<3_D, T> detector_location(int projection) const override final {
        return transform_location_(image_center_() +
                                       (relative_detector_distance_(projection) *
                                        this->volume_[0] *
                                        math::standard_basis<3_D, T>(0)),
                                   projection);
    }

    std::array<math::vec<3_D, T>, 2>
    detector_tilt(int projection) const override final {
        return {apply_rotation_(math::standard_basis<3_D, T>(1), projection),
                apply_rotation_(math::standard_basis<3_D, T>(2), projection)};
    }

  private:
    inline T relative_source_distance_(int projection) const {
        return (T)2.0 - (T)0.5 * math::cos(math::pi<T> * (T)projection / this->projection_count_);
    }

    inline T relative_detector_distance_(int projection) const {
        return relative_source_distance_(projection);
    }

    inline math::vec<3_D, T> transform_location_(math::vec<3_D, T> location,
                                                 int projection) const {
        return apply_rotation_(location - image_center_(), projection) +
               image_center_();
    }

    inline math::vec<3_D, T> apply_rotation_(math::vec<3_D, T> location,
                                             int projection) const {
        static auto axis = math::standard_basis<3_D, T>(1);

        T angle_projection = (T)2.0 * math::pi<T> / this->projection_count_;

        return math::rotate(location, axis, angle_projection * projection);
    }

    inline math::vec<3_D, T> image_center_() const {
        math::vec<3_D, T> image_center = {this->volume_[0] * (T)0.5,
                                          this->volume_[1] * (T)0.5,
                                          this->volume_[2] * (T)0.5};

        return image_center;
    }
};

} // namespace geometry
} // namespace tomo
