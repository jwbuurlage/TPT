// We hide our math (vectors, intersections, ...) in this file.
// currently the intention is to forward as much as possible to glm, but maybe
// modify the syntax for our use cases
//
// FIXME: this is becoming a mess with glm types and our own

#pragma once

#include <algorithm>
#include <iostream>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>

#include "common.hpp"
#include "volume.hpp"

namespace tomo {
namespace math {

template <typename T>
constexpr T epsilon = (T)1e-8;

template <typename T>
constexpr auto pi = glm::pi<T>();

template <typename T>
constexpr T sqrt2 = (T)1.41421356237;

// matrices are given by line distributions
template <typename T = default_scalar_type>
struct matrix_element {
    int index;
    T value;
};

// Definitions for our own vector types
template <dimension D, typename T>
struct vec_type {
    // FIXME: for higher dimensional support, add own vector class for non-glm
    // sized vectors
    /* using type = my_own_vector_class<D, T>; */
};

template <typename T>
struct vec_type<1_D, T> {
    using type = T;
};

template <typename T>
struct vec_type<2_D, T> {
    using type = glm::tvec2<T>;
};

template <typename T>
struct vec_type<3_D, T> {
    using type = glm::tvec3<T>;
};

template <typename T>
struct vec_type<4_D, T> {
    using type = glm::tvec4<T>;
};

/* Generic vector type for any dimension */
template <dimension D, typename T>
using vec = typename vec_type<D, T>::type;

/* Some shorthands for specific dimension, not that usage of this may indicate
 * non-generic code(!), so use sparingly */
template <typename T>
using vec2 = vec<2_D, T>;

template <typename T>
using vec3 = vec<3_D, T>;

template <typename T>
using vec4 = vec<4_D, T>;

/* Common math operations */
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
constexpr auto sqrt(T obj) {
    return glm::sqrt(obj);
}

template <typename T>
auto normalize(T obj) {
    return glm::normalize(obj);
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
constexpr T pow(T a, int n) {
    int result = a;
    // FIXME use exponential powering
    for (int i = 1; i < n; ++i) {
        result *= a;
    }
    return result;
}

/* compute the intersection of the line p --> p2, with the line q --> q2, or
 * returning origin if there is no such point
 * TODO: return pair of boolean here better than origin as error */
template <typename T>
vec2<T> intersection(vec2<T> p, vec2<T> p2, vec2<T> q, vec2<T> q2) {
    auto r = p2 - p;
    auto s = q2 - q;

    auto r_cross_s = cross<T>(r, s);

    if (r_cross_s > epsilon<T>) {
        auto t = cross<T>(q - p, s) / r_cross_s;
        auto u = cross<T>(q - p, r) / r_cross_s;
        if (0.0 <= t && t <= 1.0 && 0.0 <= u && u <= 1.0) {
            return p + t * r;
        }
    }

    return vec2<T>(0, 0);
}

/* FIXME: we assume a more or less square box */
template <typename T>
vec2<T> box_intersection(typename vec2<T>::type p, typename vec2<T>::type p2,
                         typename vec2<T>::type box) {
    static const std::array<std::array<math::vec2<T>, 2>, 4> lines = {
        std::array<math::vec2<T>, 2>{math::vec2<T>(0, 0),
                                     math::vec2<T>(0, box.y)},
        std::array<math::vec2<T>, 2>{math::vec2<T>(0, box.y),
                                     math::vec2<T>(box.x, box.y)},
        std::array<math::vec2<T>, 2>{math::vec2<T>(box.x, 0),
                                     math::vec2<T>(box.x, box.y)},
        std::array<math::vec2<T>, 2>{math::vec2<T>(0, 0),
                                     math::vec2<T>(box.x, 0)}};

    T min_distance = std::numeric_limits<T>::max();
    math::vec2<T> best_point;
    best_point.x = -1;
    for (auto& line_segment : lines) {
        auto intersection_point =
            math::intersection<T>(p, p2, line_segment[0], line_segment[1]);
        if (intersection_point != math::vec2<T>(0, 0)) {
            auto dist = math::distance<2_D, T>(p, intersection_point);
            if (dist < min_distance) {
                best_point = intersection_point;
                min_distance = dist;
            }
        }
    }

    return best_point;
}

/* Check whether a vector `a` lies in the *open* box defined by `vol`. */
template <dimension D, typename T>
bool inside(vec<D, T> a, volume<D> vol) {
    for (int dim = 0; dim < D; ++dim) {
        if (a[dim] <= (T)0 || a[dim] >= (T)vol[dim])
            return false;
    }
    return true;
}

template <dimension D, typename T>
void interpolate(vec<D, T> a, volume<D> vol,
                 std::vector<math::matrix_element<T>>& queue) {
    // First we see what cell corner we are closest to
    vec<D, int> b = floor(a + vec<D, T>(0.5));

    // the corresponding indices for each dimension
    std::array<std::array<int, 2>, D> cells_indices;
    for (int i = 0; i < D; ++i) {
        cells_indices[i] = {b[i] - 1, b[i]};
    }

    /* This cartesian product uses bit shifts to compute it efficiently, we can
     * do this because each set only has two elements. Note: requires D < 32,
     * which we should never exceed, ever. */
    auto cartesian_product = [](std::array<std::array<int, 2>, D> indices) {
        std::array<vec<D, int>, pow(2, D)> result;
        for (int k = 0; k < (1 << D); ++k) {
            for (int i = 0; i < D; ++i) {
                result[k][i] = indices[i][(k & (1 << i)) != 0];
            }
        }
        return result;
    };
    auto cells = cartesian_product(cells_indices);

    // Next we do a general interpolation
    for (auto cell : cells) {
        vec<D, T> cell_center = vec<D, T>(cell) + vec<D, T>((T)0.5);
        if (inside<D, T>(cell_center, vol)) {
            int index = vol.index_by_vector(cell);
            T value = (T)1.0 - ((distance(a, cell_center)) / sqrt<T>(D));
            queue.push_back({index, value});
        }
    }
}

} // namespace math
} // namespace tomo
