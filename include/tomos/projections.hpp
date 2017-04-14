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
    projections(const geometry::base<D, T>& geometry)
        : projection_count_(geometry.projection_count()),
          shape_(geometry.detector_shape()),
          projection_pixels_(geometry.detector_pixel_count()) {
        data_.resize(geometry.lines());
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

    /** Clear the projection stack. Sets each measurement to zero. */
    void clear() { std::fill(data_.begin(), data_.end(), 0); }

    /** Get a projection as an image */
    image<D - 1, T> get_projection(int projection) const {
        auto img = image<D - 1, T>(volume<D - 1, T>(shape_));
        for (int i = 0; i < projection_pixels_; ++i) {
            img[i] = data_[projection_pixels_ * projection + i];
        }
        return img;
    }

    /** Set a projection using an image */
    void set_projection(int projection, const image<D - 1, T>& img) {
        for (int i = 0; i < projection_pixels_; ++i) {
            data_[projection_pixels_ * projection + i] = img[i];
        }
    }

  private:
    int projection_count_;
    math::vec<D - 1, int> shape_;
    int projection_pixels_;
    std::vector<T> data_;
};

} // namespace tomo
