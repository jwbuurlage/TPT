#pragma once

#include <limits>

#include "../geometry.hpp"
#include "../math.hpp"
#include "../volume.hpp"
#include "trajectory.hpp"

namespace tomo {
namespace geometry {

/**
 * A laminographic geometry
 *
 * \tparam D the dimension of the volume.
 * \tparam T the scalar type to use
 */
template <typename T>
class laminography : public trajectory<3_D, T> {
  public:
    laminography(volume<3_D, T> volume, int projection_count,
                 math::vec<2_D, T> detector_size,
                 math::vec<2_D, int> detector_shape,
                 math::vec<3_D, T> source_position,
                 math::vec<3_D, T> detector_position, T source_radius,
                 T detector_radius)
        : trajectory<3_D, T>(volume, projection_count, detector_size,
                             detector_shape),
          source_position_(source_position),
          detector_position_(detector_position), source_radius_(source_radius),
          detector_radius_(detector_radius) {}

    laminography(volume<3_D, T> volume, int projection_count,
                 math::vec<2_D, T> detector_size,
                 math::vec<2_D, int> detector_shape,
                 T source_to_center = (T)1.0, T detector_to_center = (T)1.0,
                 T source_radius = (T)1.0, T detector_radius = (T)1.0)
        : laminography(volume, projection_count, detector_size, detector_shape,
                       math::volume_center(volume) -
                           source_to_center * math::standard_basis<3_D, T>(2),
                       math::volume_center(volume) +
                           detector_to_center * math::standard_basis<3_D, T>(2),
                       source_radius, detector_radius) {}

    math::vec<3_D, T> source_location(int projection) const override final {
        return source_position_ +
               apply_rotation_(source_radius_ * math::standard_basis<3_D, T>(0),
                               projection);
    }

    math::vec<3_D, T> detector_location(int projection) const override final {
        return detector_position_ +
               apply_rotation_(-detector_radius_ *
                                   math::standard_basis<3_D, T>(0),
                               projection);
    }

    std::array<math::vec<3_D, T>, 2>
    detector_tilt(int projection) const override final {
        auto a = math::standard_basis<3_D, T>(2);
        auto b = math::normalize(detector_location(projection) -
                                 source_location(projection));
        auto n = math::cross<T>(a, b);
        auto n_norm = math::norm<3_D, T>(n);
        auto theta = -math::asin<T>(n_norm);
        auto n_hat = n / n_norm;

        return {math::rotate(math::standard_basis<3_D, T>(0), n_hat, theta) *
                    this->detector_size_[0] * (T)0.5,
                math::rotate(math::standard_basis<3_D, T>(1), n_hat, theta) *
                    this->detector_size_[1] * (T)0.5};
    }

    T& source_radius() { return source_radius_; }
    T& detector_radius() { return detector_radius_; }

  private:
    math::vec<3_D, T> source_position_;
    math::vec<3_D, T> detector_position_;

    T source_radius_;
    T detector_radius_;

    inline math::vec<3_D, T> transform_location_(math::vec<3_D, T> location,
                                                 int projection) const {
        return apply_rotation_(location - image_center_(), projection) +
               image_center_();
    }

    inline math::vec<3_D, T> apply_rotation_(math::vec<3_D, T> location,
                                             int projection) const {
        static auto axis = math::standard_basis<3_D, T>(2);

        T angle_projection = (T)2.0 * math::pi<T> / this->projection_count_;

        return math::rotate(location, axis, angle_projection * projection);
    }

    inline math::vec<3_D, T> image_center_() const {
        math::vec<3_D, T> image_center =
            this->volume_.origin() + (T)0.5 * this->volume_.physical_lengths();

        return image_center;
    }
};

} // namespace geometry
} // namespace tomo
