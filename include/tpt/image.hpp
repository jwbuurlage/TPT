#pragma once

#include <algorithm>
#include <stddef.h>
#include <type_traits>
#include <vector>

#include "common.hpp"
#include "volume.hpp"

namespace tpt {

/**
 * A discrete representation of the (reconstructed) object.
 *
 * An image object holds the values for the individual voxels within a volume.
 *
 * \tparam D the dimension of the volume (and thus the reconstruction problem).
 * \tparam T the scalar type to use
 */
template <dimension D, typename T = default_scalar_type>
class image {
  public:
    using value_type = T;

    /** Construct a default-initialized image for a given volume. */
    image(volume<D, T> v) : v_(v), data_(v.cells()) {}

    /** Construct a image for a given volume with a constant value. */
    image(volume<D, T> v, T value) : v_(v), data_(v.cells(), value) {}

    /** Obtain the index of an image voxel within the volume. */
    size_t index(math::vec<D, int> xs) const { return v_.index(xs); }

    /**
     * Obtain the value for the an image voxel.
     *
     * \param index the global index of the voxel
     */
    T& operator[](size_t index) { return data_[index]; }
    const T& operator[](size_t index) const { return data_[index]; }

    /**
     * Obtain the value for the an image voxel.
     *
     * \param xs an array containing the indices of the voxel in each dimension.
     */
    T& operator()(math::vec<D, int> xs) { return (*this)[index(xs)]; }
    const T& operator()(math::vec<D, int> xs) const {
        return (*this)[index(xs)];
    }

    /**
     * Obtain the size of the image for a given dimension
     * \param i the dimension to query
     */
    int size(int i) const { return v_.voxels()[i]; }
    auto size() const { return data_.size(); }

    /** Obtain a reference to the underlying image data. */
    std::vector<T>& mutable_data() { return data_; }
    const std::vector<T>& data() const { return data_; }

    /** Obtain the volume. */
    volume<D, T> get_volume() const { return v_; }

    /** Clear the image. Fills each voxel with zero. */
    void clear() { std::fill(data_.begin(), data_.end(), 0); }

    /** Obtain an iterator to the first voxel of the image. */
    auto begin() { return data_.begin(); }

    /** Obtain an iterator that points beyond the last first voxel. */
    auto end() { return data_.end(); }

    /** Return the image shape */
    auto shape() { return v_.voxels(); }

  private:
    volume<D, T> v_;
    std::vector<T> data_;
};

} // namespace tpt
