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
    using super = geometry::trajectory<3_D, T>;
    using super::super;

    math::vec<3_D, T> source_location(int step) const override final {
        return transform_location_(
            -(T)1.0 * this->volume_[0] * math::standard_basis<3_D, T>(0), step);
    }

    math::vec<3_D, T> detector_location(int step) const override final {
        return transform_location_(
            (T)2.0 * this->volume_[0] * math::standard_basis<3_D, T>(0), step);
    }

    std::array<math::vec<3_D, T>, 2>
    detector_tilt(int step) const override final {
        return {transform_location_(math::standard_basis<3_D, T>(1), step),
                transform_location_(math::standard_basis<3_D, T>(2), step)};
    }

  private:
    inline math::vec<3_D, T> transform_location_(math::vec<3_D, T> location,
                                                 int step) const {
        static std::array<math::vec<3_D, T>, 2> axes_ = {
            math::standard_basis<3_D, T>(1), math::standard_basis<3_D, T>(2)};

        auto current_axis = axes_[step / (this->steps_ / 2)];
        auto current_step = step % (this->steps_ / 2);
        T angle_step = (T)current_step / (math::pi<T> * this->steps_);

        return math::rotate(location, current_axis, angle_step * current_step);
    }
};

} // namespace geometry
} // namespace tomo
