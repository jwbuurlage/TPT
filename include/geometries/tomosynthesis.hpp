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
    tomosynthesis(volume<3_D, T> volume, int projection_count,
                  math::vec<2_D, T> detector_size,
                  math::vec<2_D, int> detector_shape,
                  T relative_source_distance = (T)1.0,
                  T relative_detector_distance = (T)1.0,
                  T source_arc = (T)0.5 * math::pi<T>)
        : trajectory<3_D, T>(volume, projection_count, detector_size, detector_shape),
          relative_source_distance_(relative_source_distance),
          relative_detector_distance_(relative_detector_distance),
          source_arc_(source_arc) {}

    math::vec<3_D, T> source_location(int projection) const override final {
        auto pivot = image_center_() +
                     relative_source_distance_ * this->volume_[1] *
                         math::standard_basis<3_D, T>(1);

        return apply_rotation_(pivot, projection);
    }

    math::vec<3_D, T> detector_location(int /* projection */) const override final {
        auto pivot = image_center_() -
                     relative_detector_distance_ * this->volume_[1] *
                         math::standard_basis<3_D, T>(1);

        return pivot;
    }

    std::array<math::vec<3_D, T>, 2>
    detector_tilt(int /* projection */) const override final {
        return {math::standard_basis<3_D, T>(0),
                math::standard_basis<3_D, T>(2)};
    }

    T& relative_source_distance() { return relative_source_distance_; }
    T& relative_detector_distance() { return relative_detector_distance_; }
    T& source_arc() { return source_arc_; }

  private:
    T relative_source_distance_;
    T relative_detector_distance_;
    T source_arc_;

    inline math::vec<3_D, T> apply_rotation_(math::vec<3_D, T> location,
                                             int projection) const {
        static auto axis = math::standard_basis<3_D, T>(0);
        T angle_projection = source_arc_ / this->projection_count_;
        return math::rotate(location - image_center_(), axis,
                            -source_arc_ * (T)0.5 + angle_projection * projection) +
               image_center_();
    }

    inline math::vec<3_D, T> image_center_() const {
        math::vec<3_D, T> image_center =
            this->volume_.origin() + (T)0.5 * this->volume_.physical_lengths();

        return image_center;
    }
};

} // namespace geometry
} // namespace tomo
