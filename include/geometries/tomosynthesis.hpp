#pragma once

#include <limits>

#include "../geometry.hpp"
#include "../math.hpp"
#include "../volume.hpp"
#include "trajectory.hpp"

namespace tomo {
namespace geometry {

/**
 * A tomosynthesis geometry
 *
 * \tparam D the dimension of the volume.
 * \tparam T the scalar type to use
 */
template <typename T>
class tomosynthesis : public trajectory<3_D, T> {
  public:
    /** Construct the geometry with a given number of lines. */
    tomosynthesis(volume<3_D> volume, int steps, T detector_spacing = (T)1.0,
                  math::vec<2_D, int> detector_size = math::vec<2_D, int>{1},
                  T relative_source_distance = (T)1.0,
                  T relative_detector_distance = (T)1.0,
                  T source_radius = (T)1.0, T source_arc = (T)0.5 * math::pi<T>)
        : trajectory<3_D, T>(volume, steps, detector_spacing, detector_size),
          relative_source_distance_(relative_source_distance),
          relative_detector_distance_(relative_detector_distance),
          source_radius_(source_radius), source_arc_(source_arc) {}

    math::vec<3_D, T> source_location(int step) const override final {
        auto pivot = image_center_() +
                     relative_source_distance_ * this->volume_.y() *
                         math::standard_basis<3_D, T>(1);

        return apply_rotation_(pivot, step);
    }

    math::vec<3_D, T> detector_location(int step) const override final {
        auto pivot = image_center_() -
                     relative_detector_distance_ * this->volume_.y() *
                         math::standard_basis<3_D, T>(1);

        return pivot;
    }

    std::array<math::vec<3_D, T>, 2>
    detector_tilt(int step) const override final {
        return {math::standard_basis<3_D, T>(0),
                math::standard_basis<3_D, T>(2)};
    }

  private:
    T relative_source_distance_;
    T relative_detector_distance_;
    T source_radius_;
    T source_arc_;

    inline math::vec<3_D, T> apply_rotation_(math::vec<3_D, T> location,
                                             int step) const {
        static auto axis = math::standard_basis<3_D, T>(2);
        T angle_step = source_arc_ / this->steps_;
        return math::rotate(location, axis, -source_arc_ * (T)0.5 + angle_step * step);
    }

    inline math::vec<3_D, T> image_center_() const {
        math::vec<3_D, T> image_center = {this->volume_.x() * (T)0.5,
                                          this->volume_.y() * (T)0.5,
                                          this->volume_.z() * (T)0.5};

        return image_center;
    }
};

} // namespace geometry
} // namespace tomo
