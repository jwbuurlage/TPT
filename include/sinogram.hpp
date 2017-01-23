#pragma once

#include <vector>

#include "common.hpp"

namespace tomo {

/**
 * FIXME
 * - I want to redo this.
 * - [ ] Definition of a sinogram
 * - [ ] Stack of projections
 * - [ ] What about 2d vs 3d vs Dd?
 */

/**
 * An object that represents the measurements for a given geometry (and
 * projector).
 *
 * \tparam D the dimension of the volume (and thus the reconstruction problem).
 * \tparam T the scalar type to use
 * \tparam Geometry the geometry corresponding to the measurements
 * \tparam Projector the projector corresponding to the measurement. This
 * particular template parameter can be used by reconstruction algorithms to
 * compare the sinogram using the appropriate ray-integration method.
 */
template <dimension D, typename T, class Geometry, class Projector>
class sinogram {
  public:
    /** The type of the projector used for this sinogram. */
    using projector_type = Projector;
    using value_type = T;

    /** Construct a default-initialized sinogram for a geometry. */
    sinogram(const Geometry& geometry) : geometry_(geometry) {
        data_.resize(geometry_.lines());
    }

    /**
     * Obtain a reference to the i-th measurement, corresponding to the i-th
     * line.
     * \param i the index of the line
     */
    T& operator[](int i) { return data_[i]; }
    const T& operator[](int i) const { return data_[i]; }

    /**
     * Obtain a reference to the underlying sinogram data.
     */
    const std::vector<T>& data() const { return data_; }
    std::vector<T>& mutable_data() { return data_; }

    /** Obtain a reference to the geometry. */
    const Geometry& geometry() const { return geometry_; }

    /** Clear the sinogram. Sets each measurement to zero. */
    void clear() { std::fill(data_.begin(), data_.end(), 0); }

    auto get_volume() const { return geometry().get_volume(); }

    image<D, T> as_image() const {
        image<2_D, T> img(
            volume<2_D>(geometry_.groups()[0], geometry_.groups()[1]));
        img.mutable_data() = data_;
        return img;
    }

  private:
    const Geometry& geometry_;
    std::vector<T> data_;
};

} // namespace tomo
