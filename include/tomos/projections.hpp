#pragma once

#include <algorithm>
#include <numeric>
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
        : geometry_(geometry), data_(geometry.lines()) {
        compute_offsets_();
    }

    /** Construct projections for a geometry with constant value. */
    projections(const geometry::base<D, T>& geometry, T value)
        : geometry_(geometry), data_(geometry.lines(), value) {
        compute_offsets_();
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

    ///** Get a projection as an image */
    image<D - 1, T> get_projection(int idx) const {
        auto offset = offsets_[idx];
        auto shape = geometry_.projection_shape(idx);
        auto img = image<D - 1, T>(volume<D - 1, T>(shape));
        for (int i = 0; i < math::reduce<D - 1>(shape); ++i) {
            img[i] = data_[offset + i];
        }
        return img;
    }

    /** Set a projection using an image */
    void set_projection(int idx, const image<D - 1, T>& img) {
        auto offset = offsets_[idx];
        auto pixels = math::reduce<D - 1>(geometry_.projection_shape(idx));
        for (int i = 0; i < pixels; ++i) {
            data_[offset + i] = img[i];
        }
    }

    int offset(int idx) const { return offsets_[idx]; }

  private:
    void compute_offsets_() {
        offsets_.resize(geometry_.projection_count());
        offsets_[0] = 0;
        std::iota(offsets_.begin() + 1, offsets_.end(), 0);
        std::transform(offsets_.begin() + 1, offsets_.end(),
                       offsets_.begin() + 1, [&](int i) {
                           return math::reduce<D - 1>(
                               geometry_.projection_shape(i));
                       });
        std::partial_sum(offsets_.begin() + 1, offsets_.end(),
                         offsets_.begin() + 1);
    }

    std::vector<int> offsets_;
    const geometry::base<D, T>& geometry_;
    int projection_count_;
    std::vector<T> data_;
};

} // namespace tomo
