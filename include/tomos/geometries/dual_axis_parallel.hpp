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
    dual_axis_parallel(volume<3_D, T> volume, int projection_count_per_axis,
                       math::vec<2_D, T> detector_size,
                       math::vec<2_D, int> detector_shape)
        : trajectory<3_D, T>(volume, 2 * projection_count_per_axis,
                             detector_size, detector_shape),
          actual_detector_size_(detector_size) {
        static_assert(false, "(dual-axis) parallel support should be "
                             "implemented with a 'parallel trajectory' base");
        projection_count_per_axis_ = math::reduce<2_D>(detector_size);
        assert(projection_count_per_axis_ > 0);
    }

    math::vec<3_D, T> source_location(int projection) const override final {
        return transform_location_(
            math::volume_center(this->volume_) -
                ((T) this->volume_[0] * math::standard_basis<3_D, T>(0)),
            projection_actual_(projection));
    }

    math::vec<3_D, T> detector_location(int projection) const override final {
        return transform_location_(
                   math::volume_center(this->volume_) +
                       ((T) this->volume_[0] * math::standard_basis<3_D, T>(0)),
                   projection_actual_(projection)) +
               this->detector_offset_(projection,
                                      projection % projection_count_per_axis_,
                                      actual_detector_size_);
    }

    std::array<math::vec<3_D, T>, 2>
    detector_tilt(int projection) const override final {
        return {apply_rotation_(math::standard_basis<3_D, T>(1),
                                projection_actual_(projection)),
                apply_rotation_(math::standard_basis<3_D, T>(2),
                                projection_actual_(projection))};
    }

  private:
    inline math::vec<3_D, T> transform_location_(math::vec<3_D, T> location,
                                                 int projection) const {
        return apply_rotation_(location - math::volume_center(this->volume_),
                               projection) +
               math::volume_center(this->volume_);
    }

    inline math::vec<3_D, T> apply_rotation_(math::vec<3_D, T> location,
                                             int projection) const {
        static auto axes = std::array<math::vec<3_D, T>, 2>{
            math::standard_basis<3_D, T>(1), math::standard_basis<3_D, T>(2)};
        auto axis = axes[projection / projection_count_per_axis_];

        T angle_projection = math::pi<T> / projection_count_per_axis_;

        return math::rotate(location, axis,
                            angle_projection *
                                (projection % projection_count_per_axis_));
    }

    int projection_actual_(int projection) const {
        return projection / projection_count_per_axis_;
    }

    int projection_count_per_axis_;
    math::vec<2_D, int> actual_detector_size_;
};

} // namespace geometry
} // namespace tomo
