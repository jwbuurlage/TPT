#pragma once

#include <array>
#include <cassert>
#include <vector>

#include "../geometry.hpp"
#include "../math.hpp"

namespace tomo {
namespace geometry {

/**
 * A geometry defined by trajectories of the source and the detector (array)
 *
 * The detector and the source move in projections along some trajectory, which
 * is described by concrete subclasses.
 *
 * \tparam D the dimension of the volume.
 * \tparam T the scalar type to use
 */
template <dimension D, typename T>
class trajectory : public base<D, T> {
  public:
    /** Construct the geometry with a given number of lines. */
    trajectory(volume<D, T> volume, int projection_count,
               math::vec<D - 1, T> detector_size,
               math::vec<D - 1, int> detector_shape)
        : base<D, T>(projection_count), volume_(volume),
          detector_size_(detector_size), detector_shape_(detector_shape) {
        this->compute_lines_();
    }

    /**
     * The tilt of the detector in projection `projection`. The tilt is given as
     * the *normalized* principal axes of the hyperplane, i.e. by (D - 1)
     * D-dimensional vectors.
     */
    virtual std::array<math::vec<D, T>, D - 1>
    detector_tilt(int projection) const = 0;

    auto detector_size() const { return detector_size_; }

    virtual math::vec<D, T> detector_location(int i) const = 0;

    math::vec<D, T> detector_corner(int i) const override {
        math::vec<D, T> result = detector_location(i);
        for (int d = 0; d < D - 1; ++d) {
            result -= (T)0.5 * (detector_size_[d] * detector_tilt(i)[d]);
        }
        return result;
    }

    math::vec<D - 1, int> projection_shape(int) const override {
        return detector_shape_;
    }

    projection<D, T> get_projection(int idx) const {
        return {this->source_location(idx),
                this->detector_location(idx),
                detector_size_,
                detector_tilt(idx),
                this->detector_shape_,
                false};
    }

    std::array<math::vec<D, T>, D - 1> projection_delta(int i) const override {
        auto axes = detector_tilt(i);
        std::array<math::vec<D, T>, D - 1> result;

        for (int d = 0; d < D - 1; ++d) {
            result[d] =
                (detector_size_[d] / this->detector_shape_[d]) * axes[d];
        }

        return result;
    }

  protected:
    virtual ~trajectory() = default;

    volume<D, T> volume_;

    math::vec<D - 1, T> detector_size_;
    math::vec<D - 1, int> detector_shape_;
};

} // namespace geometry
} // namespace tomo
