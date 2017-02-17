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
class cone_beam : public trajectory<3_D, T> {
  public:
    /** Construct the geometry with a given number of lines. */
    cone_beam(volume<3_D, T> volume, int projection_count,
              math::vec<2_D, T> detector_size,
              math::vec<2_D, int> detector_shape,
              math::vec<3_D, T> source_position,
              math::vec<3_D, T> detector_position,
              std::array<math::vec<3_D, T>, 2> detector_tilt,
              std::vector<T> angles = {})
        : trajectory<3_D, T>(volume, projection_count, detector_size,
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

    math::vec<3_D, T> source_location(int projection) const override final {
        return transform_location_(source_position_, projection);
    }

    math::vec<3_D, T> detector_location(int projection) const override final {
        return transform_location_(detector_position_, projection);
    }

    std::array<math::vec<3_D, T>, 2>
    detector_tilt(int projection) const override final {
        return {apply_rotation_(detector_tilt_[0], projection),
                apply_rotation_(detector_tilt_[1], projection)};
    }

    void reduce_projections(int scaled_projection_count) override {
        auto step = (T)(this->projection_count_) / scaled_projection_count;
        for (int i = 0; i < scaled_projection_count; ++i) {
            int j = (int)(i * step);
            angles_[i] = angles_[j];
        }
        angles_.resize(scaled_projection_count);
        this->projection_count_ = scaled_projection_count;
        this->line_count_ =
            this->projection_count_ * this->detector_pixel_count_;
    }

  private:
    math::vec<3_D, T> source_position_;
    math::vec<3_D, T> detector_position_;
    std::array<math::vec<3_D, T>, 2> detector_tilt_;
    std::vector<T> angles_;

    inline math::vec<3_D, T> transform_location_(math::vec<3_D, T> location,
                                                 int projection) const {
        return apply_rotation_(location - math::volume_center(this->volume_),
                               projection) +
               math::volume_center(this->volume_);
    }

    inline math::vec<3_D, T> apply_rotation_(math::vec<3_D, T> location,
                                             int projection) const {
        // here we rotate around z
        return math::rotate(location, math::standard_basis<3_D, T>(2),
                            angles_[projection]);
    }
};

} // namespace geometry
} // namespace tomo
