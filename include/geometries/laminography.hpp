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
    /** Construct the geometry with a given number of lines. */
    laminography(volume<3_D> volume, int steps, T detector_spacing = (T)1.0,
                 math::vec<2_D, int> detector_size = math::vec<2_D, int>{1},
                 T relative_source_distance = (T)1.0,
                 T relative_detector_distance = (T)1.0,
                 T source_radius = (T)1.0, T detector_radius = (T)1.0)
        : trajectory<3_D, T>(volume, steps, detector_spacing, detector_size),
          relative_source_distance_(relative_source_distance),
          relative_detector_distance_(relative_detector_distance),
          source_radius_(source_radius), detector_radius_(detector_radius) {}

    math::vec<3_D, T> source_location(int step) const override final {
        auto pivot = image_center_() -
                     relative_source_distance_ * this->volume_.y() *
                         math::standard_basis<3_D, T>(1);

        return pivot +
               apply_rotation_(source_radius_ * math::standard_basis<3_D, T>(0),
                               step);
    }

    math::vec<3_D, T> detector_location(int step) const override final {
        auto pivot = image_center_() +
                     relative_detector_distance_ * this->volume_.y() *
                         math::standard_basis<3_D, T>(1);

        return pivot +
               apply_rotation_(
                   -detector_radius_ * math::standard_basis<3_D, T>(0), step);
    }

    std::array<math::vec<3_D, T>, 2>
    detector_tilt(int step) const override final {
        auto a = math::standard_basis<3_D, T>(1);
        auto b =
            math::normalize(detector_location(step) - source_location(step));
        auto n = math::cross<T>(a, b);
        auto n_norm = math::norm<3_D, T>(n);
        auto theta = math::asin<T>(n_norm);
        auto n_hat = n / n_norm;

        return {math::rotate(math::standard_basis<3_D, T>(0), n_hat, theta),
                math::rotate(math::standard_basis<3_D, T>(2), n_hat, theta)};
    }

    T& relative_source_distance() { return relative_source_distance_; }
    T& relative_detector_distance() { return relative_detector_distance_; }
    T& source_radius() { return source_radius_; }
    T& detector_radius() { return detector_radius_; }

  private:
    T relative_source_distance_;
    T relative_detector_distance_;
    T source_radius_;
    T detector_radius_;

    inline math::vec<3_D, T> transform_location_(math::vec<3_D, T> location,
                                                 int step) const {
        return apply_rotation_(location - image_center_(), step) +
               image_center_();
    }

    inline math::vec<3_D, T> apply_rotation_(math::vec<3_D, T> location,
                                             int step) const {
        static auto axis = math::standard_basis<3_D, T>(1);

        T angle_step = (T)2.0 * math::pi<T> / this->steps_;

        return math::rotate(location, axis, angle_step * step);
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
