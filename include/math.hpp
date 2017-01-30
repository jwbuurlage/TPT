/**
 * \file We hide our math (vectors, intersections, ...) in this file.
 * currently the intention is to forward as much as possible to glm, but maybe
 * modify the syntax for our use cases, and to remain flexible w.r.t.
 * alternative vector libraries.
 */

#pragma once

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

#include <experimental/optional>

template <typename T>
using optional = std::experimental::optional<T>;

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtx/rotate_vector.hpp>

// this allows for `glm::to_string(vecD)`, which eases debugging
#include <glm/gtx/string_cast.hpp>

#include "common.hpp"
#include "math/vector.hpp"
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

/** Compute the absolute value of an object of type T. */
template <typename T>
auto abs(T obj) {
    return glm::abs(obj);
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

/**
 * A ray is a line segment from source to a detector
 *
 * \tparam D the dimension of the volume.
 * \tparam T the scalar type to use
 */
template <dimension D, typename T = default_scalar_type>
struct ray {
    /** Construct a line using defaults. */
    ray() = default;

    /** Construct a line from two vectors. */
    ray(vec<D, T> source_in, vec<D, T> detector_in)
        : source(source_in), detector(detector_in) {}

    vec<D, T> source;   //> the position of the source
    vec<D, T> detector; //> the position of the detector pixel
};

/**
 * A line has an origin, and a delta -- and is used for numerical integration
 *
 * \tparam D the dimension of the volume.
 * \tparam T the scalar type to use
 */
template <dimension D, typename T = default_scalar_type>
struct line {
    /** Construct a line using defaults. */
    line() = default;

    /** Construct a line from two vectors. */
    line(vec<D, T> origin_in, vec<D, T> delta_in)
        : origin(origin_in), delta(delta_in) {}

    vec<D, T> origin; //> first intersection of line with a volume
    vec<D, T> delta;  //> the direction vector
};

/** Return the `axis`-th vector of the standard basis */
template <dimension D, typename T>
auto standard_basis(int axis) {
    auto result = vec<D, T>{};
    result[axis] = (T)1;
    return result;
}

/** Return a vector with the sign of each component as 1 or -1 */
template <dimension D, typename T>
vec<D, int> sign(vec<D, T> rhs) {
    auto result = vec<D, int>{};
    for (int d = 0; d < D; ++d) {
        if (approx_equal(rhs[d], (T)0)) {
            result[d] = 0;
        } else {
            result[d] = rhs[d] > 0 ? 1 : -1;
        }
    }
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

/** Compute the exponential of an object of type T. */
template <typename T>
auto exp(T obj) {
    return glm::exp(obj);
}

/** Floor an object of type T. */
template <typename T>
auto floor(T obj) {
    return glm::floor(obj);
}

/** Ceil an object of type T. */
template <typename T>
auto ceil(T obj) {
    return glm::ceil(obj);
}

/** Round an object of type T to the nearest integer. */
template <typename T>
auto round(T obj) {
    return glm::round(obj);
}

/** Rotate a 3D vector */
template <typename T>
auto rotate(vec<3_D, T> v, vec<3_D, T> normal, T angle) {
    return glm::rotate(v, angle, normal);
}

/** Multiply the elements of a vector together */
template <dimension D, typename T>
auto product(vec<D, T> vec) {
    auto result = (T)1;
    for (int d = 0; d < D; ++d) {
        result *= vec[d];
    }
    return result;
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

/** Return the value of the maximum component. */
template <dimension D, typename T>
T max_element(vec<D, T> v) {
    auto max = std::numeric_limits<T>::min();
    for (int d = 0; d < D; ++d) {
        if (v[d] > max) {
            max = v[d];
        }
    }
    return max;
}

/** Return the index of the maximum component. */
template <dimension D, typename T>
int max_index(vec<D, T> x) {
    int max_index = -1;
    auto max = std::numeric_limits<T>::min();
    for (int d = 0; d < D; ++d) {
        if (x[d] > max) {
            max = x[d];
            max_index = d;
        }
    }
    return max_index;
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

/** Checks whether a world vector lies in the box defined by a volume. This is
 * open by one side, the upper side in each axis. */
template <dimension D, typename T>
bool inside(vec<D, T> a, volume<D> vol) {
    a -= vol.origin();
    for (int dim = 0; dim < D; ++dim) {
        if (a[dim] < -epsilon<T> || a[dim] > (T)vol[dim] + epsilon<T>)
            return false;
    }
    return true;
}

/** Checks whether a world vector lies in the box defined by a volume. This is
 * open by one side, the upper side in each axis. */
template <dimension D, typename T>
bool inside_margin(vec<D, T> a, volume<D> vol, T margin) {
    a -= vol.origin();
    for (int dim = 0; dim < D; ++dim) {
        if (a[dim] < -(margin + epsilon<T>) ||
            a[dim] > (T)vol[dim] + margin + epsilon<T>)
            return false;
    }
    return true;
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

// intersect the line from 'a' to 'b' with the AABB cornered in the origin of
// size 'sides'
template <dimension D, typename T>
optional<std::pair<vec<D, T>, vec<D, T>>>
aabb_intersection(vec<D, T> a, vec<D, T> b, vec<D, T> sides,
                  vec<D, T> origin = vec<D, T>{0}) {
    T t_min = std::numeric_limits<T>::min();
    T t_max = std::numeric_limits<T>::max();

    a -= origin;
    b -= origin;

    for (int d = 0; d < D; ++d) {
        if (!approx_equal(a[d], b[d])) {
            auto delta = b[d] - a[d];
            auto t1 = -a[d] / delta;
            auto t2 = (sides[d] - a[d]) / delta;
            t_min = std::max(t_min, (t1 < t2) ? t1 : t2);
            t_max = std::min(t_max, (t1 > t2) ? t1 : t2);
        } else {
            if (a[d] < 0 || a[d] >= sides[d]) {
                return optional<std::pair<vec<D, T>, vec<D, T>>>();
            }
        }
    }

    if (t_max <= t_min) {
        // no intersection
        return optional<std::pair<vec<D, T>, vec<D, T>>>();
    }

    a += origin;
    b += origin;

    return optional<std::pair<vec<D, T>, vec<D, T>>>{
        std::pair<vec<D, T>, vec<D, T>>{{a + (b - a) * t_min},
                                        {a + (b - a) * t_max}}};
}

template <dimension D, typename T>
optional<line<D, T>> truncate_to_volume(ray<D, T> ray, volume<D> v) {
    // need line plane intersection, because the box is axis aligned (AABB) we
    // need to do a ray/AABB intersection.
    auto origin = aabb_intersection<D, T>(ray.source, ray.detector, v.lengths(),
                                          math::vec<D, T>(v.origin()));
    if (!origin) {
        return optional<line<D, T>>();
    }

    auto steps = round(distance(ray.source, origin.value().first));
    auto delta = normalize(ray.detector - ray.source);
    auto start = ray.source + steps * delta;

    return optional<line<D, T>>{line<D, T>{start, delta}};
}

/**
 * for a line and bounds, return a pair of vectors consisting
 * of the two intersection points with the bounds
 */
template <dimension D, typename T>
optional<std::pair<vec<D, T>, vec<D, T>>>
intersect_bounds(math::ray<D, T> l, std::array<math::vec2<int>, D> bounds) {
    vec<D, T> bounds_origin{bounds[0][0], bounds[1][0], bounds[2][0]};

    vec<D, T> bounds_sides =
        vec<D, T>{bounds[0][1], bounds[1][1], bounds[2][1]} - bounds_origin;

    auto result = math::aabb_intersection<D, T>(l.source, l.detector,
                                                bounds_sides, bounds_origin);

    if (result) {
        return std::pair<vec<D, T>, vec<D, T>>{result.value().first,
                                               result.value().second};
    }

    return result;
}

/** Reverse-interpolate a world vector to the
 * surrounding voxels. */
template <dimension D, typename T>
void interpolate(vec<D, T> a, volume<D> vol,
                 std::vector<matrix_element<T>>& queue) {
    auto relative = a - vec<D, T>(vol.origin());
    // First we see what cell corner we are closest to
    vec<D, int> b = round(relative);

    // the corresponding indices for each dimension
    std::array<std::array<int, 2>, D> cells_indices;
    for (int i = 0; i < D; ++i) {
        cells_indices[i] = {b[i] - 1, b[i]};
    }

    /* This cartesian product uses bit shifts to compute it efficiently, we
     * can
     * do this because each set only has two elements. Note: requires D <
     * 32,
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

        bool flag = false;
        // if any of the coordinates are outside, we continue
        for (int d = 0; d < D; ++d) {
            if (cell[d] < 0 || cell[d] >= vol[d]) {
                flag = true;
                break;
            }
        }
        if (flag) {
            continue;
        }

        int index = vol.index_by_vector(cell);
        auto value = math::product<D, T>(math::vec<D, T>((T)1) -
                                         math::abs(relative - cell_center));
        queue.push_back({index, value});
    }
}

template <tomo::dimension D, typename T>
math::vec<D - 1, T> restrict(math::vec<D, T> x, int skip) {
    math::vec<D - 1, T> reduced_point;
    int k = 0;
    for (int d = 0; d < D; ++d) {
        if (d == skip)
            continue;
        reduced_point[k++] = x[d];
    }
    return reduced_point;
}

template <tomo::dimension D, typename T>
math::vec<D, T> extend(math::vec<D - 1, T> x, int axis, T value) {
    math::vec<D, T> point;
    int k = 0;
    for (int d = 0; d < D; ++d) {
        if (d == axis) {
            point[d] = value;
        } else {
            point[d] = x[k++];
        }
    }
    return point;
}

} // namespace math
} // namespace tomo
