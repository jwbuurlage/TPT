#pragma once

#include <array>
#include <cassert>
#include <vector>

#include "../geometry.hpp"
#include "../math.hpp"

namespace tomo {
namespace geometry {

/**
 * Relevant data when requesting projection 'i'
 *
 * when we have this projection, we can make a subprojection by specifying:
 * - base pixel
 * - subprojection size
 */
template <dimension D, typename T>
struct projection {
    math::vec<D, T> source_location;
    math::vec<D, T> detector_location;
    math::vec<D - 1, T> detector_size;
    std::array<math::vec<D, T>, D - 1> detector_tilt;
    math::vec<D - 1, int> detector_shape;
};

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
    trajectory(volume<3_D, T> volume, int projection_count,
               math::vec<2_D, T> detector_size,
               math::vec<2_D, int> detector_shape)
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
        return detector_location(i) -
               (T)0.5 * (detector_size_[0] * detector_tilt(i)[0] +
                         detector_size_[1] * detector_tilt(i)[1]);
    }

    math::vec<D - 1, int> projection_shape(int) const override {
        return detector_shape_;
    }

    projection<D, T> get_projection(int idx) {
        return {this->source_location(idx), this->detector_location(idx),
                detector_size_, detector_tilt(idx), this->detector_shape_};
    }

    std::array<math::vec<D, T>, D - 1> projection_delta(int i) const override {
        auto axes = detector_tilt(i);
        auto relative =
            math::vec<D - 1, T>(detector_size_[0] / this->detector_shape_[0],
                                detector_size_[1] / this->detector_shape_[1]);
        return {relative[0] * axes[0], relative[1] * axes[1]};
    }

  protected:
    virtual ~trajectory() = default;

    volume<D, T> volume_;

    math::vec<2_D, T> detector_size_;
    math::vec<2_D, int> detector_shape_;
};

} // namespace geometry
} // namespace tomo
