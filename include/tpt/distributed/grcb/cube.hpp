#pragma once

#include <array>
#include <limits>

#include "../../math.hpp"
#include "datatypes.hpp"

namespace tpt::grcb {

template <typename T>
cube<T> corners(volume<3_D, T> v) {
    auto vcorners = v.corners();
    auto corners = cube<T>{};
    std::copy(vcorners.begin(), vcorners.end(), corners.begin());
    return corners;
}

template <typename T>
std::array<vec3<T>, 2> min_max_cube(cube<T> xs) {
    auto min = vec3<T>{std::numeric_limits<T>::max()};
    auto max = vec3<T>{-std::numeric_limits<T>::max()};
    for (int i = 0; i < 8; ++i) {
        for (int d = 0; d < 3; ++d) {
            if (xs[i][d] < min[d]) {
                min[d] = xs[i][d];
            }
            if (xs[i][d] > max[d]) {
                max[d] = xs[i][d];
            }
        }
    }

    return {min, max};
}

/** Complete the area of the splitting interface along the d-th axis. */
template <typename T>
T cross_section_area(cube<T> corners, size_t d) {
    auto zs = std::array<vec2<T>, 8>{};
    std::transform(corners.begin(), corners.end(), zs.begin(),
                   [=](auto x) { return math::restrict<3_D, T>(x, d); });
    auto x1 = std::numeric_limits<T>::max();
    auto x2 = -std::numeric_limits<T>::max();
    auto y1 = std::numeric_limits<T>::max();
    auto y2 = -std::numeric_limits<T>::max();
    for (auto z : zs) {
        if (z[0] < x1) {
            x1 = z[0];
        }
        if (z[0] > x2) {
            x2 = z[0];
        }
        if (z[1] < y1) {
            y1 = z[1];
        }
        if (z[1] > y2) {
            y2 = z[1];
        }
    }

    return (x2 - x1) * (y2 - y1);
}

} // namespace tpt::grcb
