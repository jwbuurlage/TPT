#pragma once

#include <array>
#include <cassert>
#include <vector>

#include "../geometry.hpp"

namespace tomo {
namespace geometry {

/**
 * A geometry defined by trajectories of the source and the detector (array)
 *
 * The detector and the source move in steps along some trajectory, which is
 * described by concrete subclasses.
 *
 * \tparam D the dimension of the volume.
 * \tparam T the scalar type to use
 */
template <dimension D, typename T>
class trajectory : public base<D, T, trajectory<D, T>> {
  public:
    /** Construct the geometry with a given number of lines. */
    trajectory(volume<D> volume, int steps, T detector_spacing = (T)1,
               math::vec<D - 1, int> detector_size = math::vec<D - 1, int>{1})
        : base<D, T, trajectory<D, T>>(steps *
                                       math::reduce<D - 1>(detector_size)),
          volume_(volume), steps_(steps), detector_size_(detector_size),
          detector_spacing_(detector_spacing),
          detector_count_(math::reduce<D - 1>(detector_size)) {}

    /**
     * Return the i-th line of this geometry.
     */
    inline line<D, T> get_line(int i) const {
        int step = i / detector_count_;
        int detector = i % detector_count_;

        auto axes = detector_tilt(step);

        // compute detector location in space
        math::vec<D, T> offset = detector_location(step);
        auto size_had = 1;
        for (int d = 0; d < D - 1; ++d) {
            auto detector_d = (detector / size_had) % detector_size_[d];
            offset += (detector_d - (detector_size_[d] - 1) * (T)0.5) *
                      detector_spacing_ * axes[d];
            size_had *= detector_size_[d];
        }

        auto source = source_location(step);
        auto delta = math::normalize(offset - source);

        return {source, delta};
    }

    /** The location of the source in step `step`. */
    virtual math::vec<D, T> source_location(int step) const = 0;

    /** The location of the detector in step `step`. */
    virtual math::vec<D, T> detector_location(int step) const = 0;

    /**
     * The tilt of the detector in step `step`. The tilt is given as the
     * principal axes of the hyperplane, i.e. by (D - 1) D-dimensional
     * vectors.
     */
    virtual std::array<math::vec<D, T>, D - 1>
    detector_tilt(int step) const = 0;

  protected:
    virtual ~trajectory() = default;

    volume<D> volume_;
    int steps_;

    math::vec<D - 1, int> detector_size_;
    T detector_spacing_;
    int detector_count_;
};

} // namespace geometry
} // namespace tomo
