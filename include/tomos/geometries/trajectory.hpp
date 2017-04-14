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
 * is
 * described by concrete subclasses.
 *
 * FIXME: maybe it is better to call this a 'focused geometry' (cone-beam
 * geometry), where parallel is an 'unfocused geometry'.
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
        : base<D, T>(projection_count, detector_shape), volume_(volume),
          detector_size_(detector_size) {}

    /**
     * Return the i-th line of this geometry.
     */
    math::ray<D, T> get_line(int i) const override final {
        int projection = i / this->detector_pixel_count_;
        int pixel = i % this->detector_pixel_count_;

        auto source = source_location(projection);
        auto target = detector_pixel_location(projection, pixel);

        return math::ray<D, T>(source, target);
    }

    /** The location of the source in projection `projection`. */
    virtual math::vec<D, T> source_location(int projection) const = 0;

    /** The location of the detector in projection `projection`. */
    virtual math::vec<D, T> detector_location(int projection) const = 0;

    math::vec<D, T> detector_pixel_location(int projection, int pixel) const {
        return detector_location(projection) +
               detector_offset_(projection, pixel);
    }

    /**
     * The tilt of the detector in projection `projection`. The tilt is given as
     * the
     * principal axes of the hyperplane, i.e. by (D - 1) D-dimensional
     * vectors.
     */
    virtual std::array<math::vec<D, T>, D - 1>
    detector_tilt(int projection) const = 0;

    auto detector_size() const { return detector_size_; }

  protected:
    virtual ~trajectory() = default;

    volume<D, T> volume_;

    math::vec<2_D, T> detector_size_;

    // compute detector location in space
    inline math::vec<D, T> detector_offset_(int projection, int pixel) const {
        math::vec<D, T> offset;
        // TODO we now recompute every pixel this location, should just have one
        // offset and loop over this projection, inefficient!
        auto axes = detector_tilt(projection);

        auto relative = math::vec<D - 1, T>(
            (T)(pixel % this->detector_shape_[0]) / this->detector_shape_[0],
            (T)(pixel / this->detector_shape_[0]) / this->detector_shape_[1]);

        offset += (relative[0] - (T)0.5) * detector_size_[0] * axes[0];
        offset += (relative[1] - (T)0.5) * detector_size_[1] * axes[1];

        return offset;
    }
};

} // namespace geometry
} // namespace tomo
