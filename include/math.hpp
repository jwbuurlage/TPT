/**
 * \file We hide our math (vectors, intersections, ...) in this file.
 * currently the intention is to forward as much as possible to glm, but maybe
 * modify the syntax for our use cases, and to remain flexible w.r.t.
 * alternative vector libraries.
 */

#pragma once

#include <algorithm>
#include <iostream>
#include <vector>

#include <experimental/optional>

template <typename T>
using optional = std::experimental::optional<T>;

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtx/rotate_vector.hpp>

#include "common.hpp"
#include "volume.hpp"

namespace tomo {

using namespace tomo::literals;

/**
 * This namespace contains all the mathematical operations used by Galactica.
 */
namespace math {

/** A constant representing a small value of type T. */
template <typename T>
constexpr T epsilon = (T)1e-5;

/** A constant representing pi with the precision of type T. */
template <typename T>
constexpr auto pi = glm::pi<T>();

/** A constant representing sqrt(2) with the precision of type T. */
template <typename T>
constexpr T sqrt2 = (T)1.41421356237;

/**
 * An element of the (implicit) precision matrix, given by an column index and
 * a value. The row is known externally.
 *
 * \tparam T the scalar type to use
 */
template <typename T = default_scalar_type>
struct matrix_element {
    /// The column index.
    int index;

    /// The value of the matrix element.
    T value;
};

/** The type to use for D-dimensional vectors. */
template <dimension D, typename T>
struct vec_type {
    // FIXME: for higher dimensional support, add own vector class for non-glm
    // sized vectors
    /* using type = my_own_vector_class<D, T>; */
};

/** The type to use for 1-dimensional vectors. */
template <typename T>
struct vec_type<1_D, T> {
    using type = std::array<T, 1>;
};

/** The type to use for 2-dimensional vectors. */
template <typename T>
struct vec_type<2_D, T> {
    using type = glm::tvec2<T>;
};

/** The type to use for 3-dimensional vectors. */
template <typename T>
struct vec_type<3_D, T> {
    using type = glm::tvec3<T>;
};

/** The type to use for 4-dimensional vectors. */
template <typename T>
struct vec_type<4_D, T> {
    using type = glm::tvec4<T>;
};

/**
 * Generic vector type for any dimension, specializing where necessary.
 * \tparam D the dimension of the vector
 * \tparam T the scalar type to use
 */
template <dimension D, typename T>
using vec = typename vec_type<D, T>::type;

/** Short-hand for 2-dimensional vectors. */
template <typename T>
using vec2 = vec<2_D, T>;

/** Short-hand for 3-dimensional vectors. */
template <typename T>
using vec3 = vec<3_D, T>;

/** Short-hand for 4-dimensional vectors. */
template <typename T>
using vec4 = vec<4_D, T>;

/**
 * A line inside a volume.
 *
 * \tparam D the dimension of the volume.
 * \tparam T the scalar type to use
 */
template <dimension D, typename T = default_scalar_type>
struct line {
    /** Construct a line using defaults. */
    line() = default;

    /** Construct a line from two vectors. */
    line(vec<D, T> origin_in, vec<D, T> delta_in, T length_in = 0)
        : origin(origin_in), delta(delta_in), length(length_in) {}

    vec<D, T> origin; //> first intersection of line with a volume
    vec<D, T> delta;  //> the direction vector
    T length = 0;
};

/** Return the `axis`-th vector of the standard basis */
template <dimension D, typename T>
auto standard_basis(int axis) {
    auto result = vec<D, T>{};
    result[axis] = (T)1;
    return result;
}

/** Compute the cosine for an object of type T. */
template <typename T>
auto cos(T obj) {
    return glm::cos(obj);
}

/** Compute the sine for an object of type T. */
template <typename T>
auto sin(T obj) {
    return glm::sin(obj);
}

/** Compute the arcsine for an object of type T. */
template <typename T>
auto asin(T obj) {
    return glm::asin(obj);
}

/** Floor an object of type T. */
template <typename T>
auto floor(T obj) {
    return glm::floor(obj);
}

/** Compute the absolute value of an object of type T. */
template <typename T>
auto abs(T obj) {
    return glm::abs(obj);
}

/** Rotate a 3D vector */
template <typename T>
auto rotate(vec<3_D, T> v, vec<3_D, T> normal, T angle) {
    return glm::rotate(v, angle, normal);
}

/** Multiply the elements of a vector together */
template <dimension D>
int reduce(vec<D, int> vec) {
    int product = 1;
    for (int d = 0; d < D; ++d) {
        product *= vec[d];
    }
    return product;
}

/** Compute the norm of a vector of type T. */
template <dimension D, typename T>
T norm(vec<D, T> vec) {
    T squared_sum = 0;
    for (int d = 0; d < D; ++d) {
        squared_sum += vec[d] * vec[d];
    }
    return sqrt(squared_sum);
}

/** Compute the square root of an object of type T. */
template <typename T>
constexpr auto sqrt(T obj) {
    return glm::sqrt(obj);
}

/** Normalize an object of type T. */
template <typename T>
auto normalize(T obj) {
    return glm::normalize(obj);
}

/** Check if two real-valued numbers are approximately equal. */
template <typename T, typename V>
bool approx_equal(T lhs, T rhs, V max_rel_diff, V max_abs_diff) {
    return abs((lhs - rhs) / rhs) <= max_rel_diff ||
           (max_abs_diff != 0 && abs(lhs - rhs) <= max_abs_diff);
}

/** Check if two real-valued numbers are approximately equal. */
template <typename T>
bool approx_equal(T lhs, T rhs) {
    return approx_equal(lhs, rhs, (T)1e-2, epsilon<T>);
}

/** Compute the distance between two vectors. */
template <dimension D, typename T>
T distance(vec<D, T> a, vec<D, T> b) {
    return glm::distance(a, b);
}

/** Compute the distance between two generic objects. */
template <typename T>
T distance(T a, T b) {
    return abs(a - b);
}

/** Compute the inner-product of two vectors. */
template <dimension D, typename T>
T dot(vec<D, T> a, vec<D, T> b) {
    return glm::dot(a, b);
}

/** Compute the cross-product of two vectors. */
template <typename T>
vec3<T> cross(vec3<T> a, vec3<T> b) {
    return glm::cross(a, b);
}

/** Compute the 'cross-product' of two 2-dimensional vectors. */
template <typename T>
T cross(typename vec<2_D, T>::type a, typename vec<2_D, T>::type b) {
    return a.x * b.y - a.y * b.x;
}

/** Compute \f$a^n\f$. */
template <typename T>
constexpr T pow(T a, int n) {
    int result = a;
    // FIXME use exponential powering
    for (int i = 1; i < n; ++i) {
        result *= a;
    }
    return result;
}

/**
 * Compute the intersection of the line \f$p \rightarrow p2\f$, with the line
 * \f$q \rightarrow q2\f$.
 *
 * \returns an optional containing a 2-dimensional vector if there is an
 * intersection, and no value otherwise.
 */
template <typename T>
optional<vec2<T>> intersection(vec2<T> p, vec2<T> p2, vec2<T> q, vec2<T> q2) {
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

    return optional<vec2<T>>();
}

/* FIXME: we assume a more or less square box */
// FIXME return optional here too?
template <typename T>
vec2<T> box_intersection(vec2<T> p, vec2<T> p2, vec2<T> box) {
    const std::array<std::array<vec2<T>, 2>, 4> lines = {
        std::array<vec2<T>, 2>{vec2<T>(0, 0), vec2<T>(0, box.y)},
        std::array<vec2<T>, 2>{vec2<T>(0, box.y), vec2<T>(box.x, box.y)},
        std::array<vec2<T>, 2>{vec2<T>(box.x, 0), vec2<T>(box.x, box.y)},
        std::array<vec2<T>, 2>{vec2<T>(0, 0), vec2<T>(box.x, 0)}};

    T min_distance = std::numeric_limits<T>::max();
    vec2<T> best_point;
    best_point.x = -1;
    for (auto& line_segment : lines) {
        auto intersection_point =
            intersection<T>(p, p2, line_segment[0], line_segment[1]);

        if (!intersection_point)
            continue;

        if (intersection_point.value() != vec2<T>(0, 0)) {
            auto dist = distance<2_D, T>(p, intersection_point.value());
            if (dist < min_distance) {
                best_point = intersection_point.value();
                min_distance = dist;
            }
        }
    }

    return best_point;
}

template <dimension D, typename T>
line<D, T> truncate_to_volume(vec<D, T> source, vec<D, T> detector, volume<D> v) {
    // need line plane intersection
    // FIXME implement
    (void)v;
    return {source, normalize(detector - source)};
}

/** Checks whether a vector lies in the *open* box defined by a volume. */
template <dimension D, typename T>
bool inside(vec<D, T> a, volume<D> vol) {
    for (int dim = 0; dim < D; ++dim) {
        if (a[dim] <= (T)0 || a[dim] >= (T)vol[dim])
            return false;
    }
    return true;
}

/** Reverse-interpolate a vector to the surrounding voxels. */
template <dimension D, typename T>
void interpolate(vec<D, T> a, volume<D> vol,
                 std::vector<matrix_element<T>>& queue) {
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
