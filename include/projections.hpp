#pragma once

#include <vector>

#include "common.hpp"

namespace tomo {

/**
 * An object that represents the measurements for a given geometry (and
 * projector).
 *
 * \tparam D the dimension of the volume (and thus the reconstruction problem).
 * \tparam T the scalar type to use
 */
template <dimension D, typename T>
class projections {
  public:
    using value_type = T;

    /** Construct default-initialized projections for a geometry. */
    projections(const geometry::base<D, T>& geometry) : geometry_(geometry) {
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
     * Obtain a reference to the underlying projection data.
     */
    const std::vector<T>& data() const { return data_; }
    std::vector<T>& mutable_data() { return data_; }

    /** Obtain a reference to the geometry. */
    const geometry::base<D, T>& geometry() const { return geometry_; }

    /** Clear the projection stack. Sets each measurement to zero. */
    void clear() { std::fill(data_.begin(), data_.end(), 0); }

    auto get_volume() const { return geometry().get_volume(); }

    image<D, T> as_image() const {
        image<2_D, T> img(
            volume<2_D, T>({geometry_.groups()[0], geometry_.groups()[1]}));
        img.mutable_data() = data_;
        return img;
    }

  private:
    const geometry::base<D, T>& geometry_;
    std::vector<T> data_;
};

} // namespace tomo
