#pragma once

#include <array>
#include <limits>
#include <vector>

#include "../common.hpp"
#include "../geometry.hpp"
#include "../math.hpp"
#include "../volume.hpp"

namespace tomo {
namespace geometry {

/**
 * Geometry defined by parallel lines with a number of views.
 *
 * \tparam D the dimension of the volume.
 * \tparam T the scalar type to use
 */
template <dimension D, typename T>
class parallel : public base<D, T> {
  public:
    /**
     * Construct the parallel geometry for a given number of angles and
     * detectors.
     *
     * \param angle_count the number of angles
     * \param volume the volume being scanned
     */
    parallel(const volume<D, T>& volume, int angle_count)
        : base<D, T>(angle_count, true), volume_(volume) {
        auto angle_step = math::pi<T> / angle_count;
        for (T angle = 0.0; angle < math::pi<T>; angle += angle_step) {
            angles_.push_back(angle);
        }

        // TODO currently we only consider cubic volumes and explicitely check
        // this here, but we can extend the definition to cuboids.
        auto is_cubic = [](tomo::volume<D, T> vol) -> bool {
            int x = vol[0];
            for (int d = 1; d < D; ++d) {
                if (vol[d] != x) {
                    return false;
                }
            }
            return true;
        };
        assert(is_cubic(volume));
        (void)is_cubic;

        source_position_ = math::volume_center(volume_);
        source_position_[0] -= (T)2.0 * volume_[0];
        detector_position_ = math::volume_center(volume_);
        detector_position_[0] += (T)2.0 * volume_[0];

        for (int d = 0; d < D - 1; ++d) {
            detector_tilt_[d] = math::standard_basis<D, T>(d + 1);
        }

        this->compute_lines_();
    }

    math::vec<D - 1, T> detector_size(int) const {
        return math::restrict<D, T>(volume_.physical_lengths(), 0);
    }

    math::vec<D - 1, int> projection_shape(int) const override {
        return math::restrict<D, int>(volume_.voxels(), D - 1);
    }

    math::vec<D, T> source_location(int i) const override final {
        return transform_location_(source_position_, i) - delta_(i);
    }

    math::vec<D, T> detector_corner(int i) const override final {
        return transform_location_(detector_position_, i) - delta_(i);
    }

    std::array<math::vec<D, T>, D - 1> detector_tilt(int i) const {
        auto result = std::array<math::vec<D, T>, D - 1>();
        for (int d = 0; d < D - 1; ++d) {
            result[d] = apply_rotation_(detector_tilt_[d], i);
        }
        return result;
    }

    std::array<math::vec<D, T>, D - 1> projection_delta(int i) const override {
        auto axes = detector_tilt(i);
        auto result = std::array<math::vec<D, T>, D - 1>();
        for (int d = 0; d < D - 1; ++d) {
            result[d] = (T)(volume_[d] / projection_shape(i)[d]) * axes[d];
        }
        return result;
    }

  private:
    volume<D, T> volume_;
    std::vector<T> angles_;

    math::vec<D, T> source_position_;
    math::vec<D, T> detector_position_;
    std::array<math::vec<D, T>, D - 1> detector_tilt_;

    inline math::vec<D, T> transform_location_(math::vec<D, T> location,
                                               int projection) const {
        return apply_rotation_(location - math::volume_center(this->volume_),
                               projection) +
               math::volume_center(this->volume_);
    }

    inline math::vec<D, T> apply_rotation_(math::vec<D, T> location,
                                           int projection) const {
        return math::rotate(location, math::standard_basis<D, T>(D - 1),
                            angles_[projection]);
    }

    inline math::vec<D, T> delta_(int i) const {
        auto result = math::vec<D, T>();
        for (int d = 0; d < D - 1; ++d) {
            result += detector_size(i)[d] * detector_tilt(i)[d];
        }
        return (T)0.5 * result;
    }
};

} // namespace geometry
} // namespace tomo
