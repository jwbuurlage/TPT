#pragma once

#include <algorithm>
#include <numeric>
#include <vector>

#include "common.hpp"

namespace tpt {

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
        : geometry_(geometry), data_(geometry.lines()) {}

    /** Construct projections for a geometry with constant value. */
    projections(const geometry::base<D, T>& geometry, T value)
        : geometry_(geometry), data_(geometry.lines(), value) {}

    /** Copy constructor */
    projections(const projections<D, T>& other)
        : geometry_(other.geometry_), data_(other.data_) {}

    /** Move constructor */
    projections(projections<D, T>&& other)
        : geometry_(other.geometry_), data_(std::move(other.data_)) {}

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

    ///** Get a projection as an image */
    image<D - 1, T> get_projection(int idx) const {
        auto offset = this->offset(idx);
        auto shape = geometry_.projection_shape(idx);
        auto img = image<D - 1, T>(volume<D - 1, T>(shape));
        for (int i = 0; i < math::reduce<D - 1>(shape); ++i) {
            img[i] = data_[offset + i];
        }
        return img;
    }

    /** Set a projection using an image */
    void set_projection(int idx, const image<D - 1, T>& img) {
        auto offset = this->offset(idx);
        auto pixels = math::reduce<D - 1>(geometry_.projection_shape(idx));
        for (int i = 0; i < pixels; ++i) {
            data_[offset + i] = img[i];
        }
    }

    int offset(int idx) const { return geometry_.offset(idx); }

    auto size() const { return geometry_.lines(); }

    auto begin() { return data_.begin(); }
    auto end() { return data_.end(); }

    /** Get the shape of the projection stack, based on the first projection size */
    auto shape() { 
        auto pixels = geometry_.projection_shape(0);
        auto result = math::vec<D, int>();
        for (auto i = 0u; i < D - 1; ++i) {
            result[i] = pixels[i];
        }
        result[D - 1] = geometry_.projection_count();
        return result;
    }

  private:
    const geometry::base<D, T>& geometry_;
    std::vector<T> data_;
};

} // namespace tpt
