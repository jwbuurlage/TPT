// We hide our math (vectors, intersections, ...) in this file.
// currently the intention is to forward as much as possible to glm, but maybe
// modify the syntax for our use cases

#pragma once

#include <algorithm>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>

#include "common.hpp"
#include "volume.hpp"

namespace tomo {
namespace math {

constexpr double epsilon = 1e-8;
constexpr auto pi = glm::pi<double>();
constexpr double sqrt2 = 1.41421356237;

// matrices are given by line distributions
template <typename T = default_scalar_type>
struct matrix_element {
    int index;
    T value;
};

// vector types
template <dimension D, typename T>
struct vec {};

template <typename T>
struct vec<2_D, T> {
    using type = glm::tvec2<T>;
};

template <typename T>
struct vec<3_D, T> {
    using type = glm::tvec3<T>;
};

template <typename T>
struct vec<4_D, T> {
    using type = glm::tvec4<T>;
};

template <typename T>
using vec2 = typename vec<2_D, T>::type;

template <typename T>
using vec3 = typename vec<3_D, T>::type;

template <typename T>
using vec4 = typename vec<4_D, T>::type;

// operations
template <typename T>
auto cos(T obj) {
    return glm::cos(obj);
}

template <typename T>
auto sin(T obj) {
    return glm::sin(obj);
}

template <typename T>
auto floor(T obj) {
    return glm::floor(obj);
}

template <typename T>
auto abs(T obj) {
    return glm::abs(obj);
}

template <typename T>
auto sqrt(T obj) {
    return glm::sqrt(obj);
}

// vector properties
template <dimension D, typename T>
T distance(typename vec<D, T>::type a, typename vec<D, T>::type b) {
    return glm::distance(a, b);
}

// vector properties
template <typename T>
T distance(T a, T b) {
    return math::abs(a - b);
}

template <dimension D, typename T>
T dot(typename vec<D, T>::type a, typename vec<D, T>::type b) {
    return glm::dot(a, b);
}

template <typename T>
vec3<T> cross(vec3<T> a, vec3<T> b) {
    return glm::cross(a, b);
}

template <typename T>
T cross(typename vec<2_D, T>::type a, typename vec<2_D, T>::type b) {
    return a.x * b.y - a.y * b.x;
}

template <typename T>
vec2<T> intersection(vec2<T> p, vec2<T> p2, vec2<T> q, vec2<T> q2) {
    auto r = p2 - p;
    auto s = q2 - q;

    auto r_cross_s = cross<T>(r, s);

    if (r_cross_s > epsilon) {
        auto t = cross<T>(q - p, s) / r_cross_s;
        auto u = cross<T>(q - p, r) / r_cross_s;
        if (0.0 <= t && t <= 1.0 && 0.0 <= u && u <= 1.0) {
            return p + t * r;
        }
    }

    return vec2<T>(0, 0);
}

template <dimension D, typename T>
bool inside(typename vec<D, T>::type a, volume<D> vol) {
    for (int dim = 0; dim < D; ++dim) {
        if (a[dim] + epsilon < 0 || a[dim] - epsilon > vol[dim])
            return false;
    }
    return true;
}

template <dimension D, typename T>
void interpolate(typename vec<D, T>::type a, volume<D> vol,
                 std::vector<math::matrix_element<T>>& queue) {
    static_assert(D == 2, "only bilinear interpolation is supported");

    // find neighbouring cells
    // we are looking for the cells with closest centers to our point
    int x = floor(a[0] + 0.5);
    int y = floor(a[1] + 0.5);

    // x-1, y    |   x,y
    //           |
    // ------- (x,y)--------
    //           |
    // x-1, y-1  |   x,y-1

    using cell_type = vec2<int>;
    cell_type cells[] = {{x - 1, y - 1}, {x, y - 1}, {x - 1, y}, {x, y}};
    for (auto cell : cells) {
        if (inside<D, double>(vec2<double>(cell) + vec2<double>(0.5), vol)) {
            int index = vol.index({cell[0], cell[1]});
            auto value =
                1.0 -
                ((distance(a, vec2<T>(cell[0] + 0.5, cell[1] + 0.5))) / sqrt2);
            queue.push_back({index, value});
        }
    }
}

} // namespace math
} // namespace tomo
