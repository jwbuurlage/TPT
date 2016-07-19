#pragma once

#include <array>

#include "common.hpp"

namespace tomo {

/**
 * An object that describes geometry of the volume that is being imaged
 */
template <dimension Dimension>
class volume {
  public:
    volume(std::array<int, Dimension> dimensions) : dimensions_(dimensions) {}

    template <typename... Ts>
    volume(Ts... dims) : dimensions_(std::array<int, Dimension>{dims...}) {}

    /**
     * Obtain the size of the first dimension
     * @return number of voxels in first dimension
     */
    int x() const { return dimensions_[0]; }

    /**
     * Obtain the size of the second dimension
     * @return number of voxels in second dimension
     * @note only valid if Dimension > 1
     */
    int y() const {
        static_assert(Dimension > 1,
                      "requesting 'y' in volume of dimension < 2");
        return dimensions_[1];
    }

    /**
     * Obtain the size of the third dimension
     * @return number of voxels in the third dimension
     * @note only valid if Dimension > 2
     */
    int z() const {
        static_assert(Dimension > 2,
                      "requesting 'z' in volume of dimension < 3");
        return dimensions_[2];
    }

    /**
     * Obtain the size of a dimension
     * @param i the index of the dimension
     * @return number of voxels in the the i-th dimension
     */
    int operator[](size_t i) const { return dimensions_[i]; }

    size_t index(std::array<int, Dimension> xs) const {
        int result = xs[0];
        int offset = dimensions_[0];
        for (int i = 1; i < Dimension; ++i) {
            result += offset * xs[i];
            offset *= dimensions_[i];
        }
        return result;
    }

    std::array<int, Dimension> dimensions() const { return dimensions_; }

  private:
    std::array<int, Dimension> dimensions_;
};

} // namespace tomo
