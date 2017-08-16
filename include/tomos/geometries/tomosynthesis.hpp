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
    tomosynthesis(volume<3_D, T> volume, int projection_count,
                  math::vec<2_D, T> detector_size,
                  math::vec<2_D, int> detector_shape,
                  math::vec<3_D, T> source_position,
                  math::vec<3_D, T> detector_position,
                  T source_arc = (T)0.7 * math::pi<T>)
        : trajectory<3_D, T>(volume, projection_count, detector_size,
                             detector_shape),
          source_position_(source_position),
          detector_position_(detector_position), source_arc_(source_arc) {}

    tomosynthesis(volume<3_D, T> volume, int projection_count,
                  math::vec<2_D, T> detector_size,
                  math::vec<2_D, int> detector_shape,
                  T source_to_center = (T)1.0, T detector_to_center = (T)1.0,
                  T source_arc = (T)0.5 * math::pi<T>)
        : tomosynthesis(volume, projection_count, detector_size, detector_shape,
                        math::volume_center(volume) +
                            source_to_center * math::standard_basis<3_D, T>(2),
                        math::volume_center(volume) -
                            detector_to_center *
                                math::standard_basis<3_D, T>(2),
                        source_arc) {}

    math::vec<3_D, T> source_location(int projection) const override final {
        return apply_rotation_(source_position_, projection);
    }

    math::vec<3_D, T>
    detector_location(int /* projection */) const override final {
        return detector_position_;
    }

    std::array<math::vec<3_D, T>, 2>
    detector_tilt(int /* projection */) const override final {
        return {math::standard_basis<3_D, T>(0),
                math::standard_basis<3_D, T>(1)};
    }

    T& source_arc() { return source_arc_; }

  private:
    math::vec<3_D, T> source_position_;
    math::vec<3_D, T> detector_position_;

    T source_arc_;

    inline math::vec<3_D, T> apply_rotation_(math::vec<3_D, T> location,
                                             int projection) const {
        static auto axis = math::standard_basis<3_D, T>(0);
        T angle_projection = source_arc_ / (this->projection_count_ - 1);
        return math::rotate(location - image_center_(), axis,
                            -source_arc_ * (T)0.5 +
                                angle_projection * projection) +
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
