#pragma once

#include <array>
#include <numeric>

#include "common.hpp"

namespace tomo {

/**
 * An object that describes geometry of the volume that is being imaged
 */
template <dimension D>
class volume {
  public:
    volume(int k) {
        for (int i = 0; i < D; ++i)
            dimensions_[i] = k;
    }

    volume(std::array<int, D> dimensions) : dimensions_(dimensions) {}

    template <typename... Ts>
    volume(Ts... dims) : dimensions_(std::array<int, D>{dims...}) {}

    /**
     * Obtain the size of the first dimension
     * @return number of voxels in first dimension
     */
    int x() const { return dimensions_[0]; }

    /**
     * Obtain the size of the second dimension
     * @return number of voxels in second dimension
     * @note only valid if D > 1
     */
    int y() const {
        static_assert(D > 1, "requesting 'y' in volume of dimension < 2");
        return dimensions_[1];
    }

    /**
     * Obtain the size of the third dimension
     * @return number of voxels in the third dimension
     * @note only valid if D > 2
     */
    int z() const {
        static_assert(D > 2, "requesting 'z' in volume of dimension < 3");
        return dimensions_[2];
    }

    /**
     * Obtain the size of a dimension
     * @param i the index of the dimension
     * @return number of voxels in the the i-th dimension
     */
    int operator[](size_t i) const { return dimensions_[i]; }

    template <typename Vector>
    int index_by_vector(Vector xs) const {
        return index_by_vector_(xs);
    }

    int index(std::array<int, D> xs) const {
        return index_by_vector_(xs);
    }

    template <typename... Ts>
    int index(Ts... xs) const {
        return index_(0, 1, xs...);
    }

    std::array<int, D> dimensions() const { return dimensions_; }

    int cells() const {
        return std::accumulate(dimensions_.begin(), dimensions_.end(), 1,
                               std::multiplies<int>());
    }

  private:
    template <typename T, typename... Ts>
    int index_(int current, int offset, T x, Ts... xs) const {
        offset *= dimensions_[D - sizeof...(xs)];
        current += offset * x;
        return index_(current, offset, xs...);
    }

    template <typename Vector>
    inline int index_by_vector_(Vector xs) const {
        int result = xs[0];
        int offset = dimensions_[0];
        for (int i = 1; i < D; ++i) {
            result += offset * xs[i];
            offset *= dimensions_[i];
        }
        return result;
    }

    int index_(int current, int offset) const { return current; }

    std::array<int, D> dimensions_;
};

} // namespace tomo
