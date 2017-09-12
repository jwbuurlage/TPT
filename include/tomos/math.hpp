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
#include <glm/gtx/rotate_vector.hpp>

// this allows for `glm::to_string(vecD)`, which eases debugging
#include <glm/gtx/string_cast.hpp>

#include "common.hpp"
#include "math/basic_operations.hpp"
#include "math/constants.hpp"
#include "math/vector.hpp"
#include "math/vector_operations.hpp"
#include "volume.hpp"

namespace tomo {

using namespace tomo::literals;

/**
 * This namespace contains all the mathematical operations used by Galactica.
 */
namespace math {

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

/** A slice is 2D section of a 3D volume, given by an orientation. */
template <typename T>
struct slice {
    /** Construct the axis-aligned slice along the d-th axis. */
    slice(int d) {
        base = vec3<T>(-1.0f);
        base[d] = 0.0f;

        b_x = (d == 0) ? vec3<T>(0.0f, 2.0f, 0.0f) : vec3<T>(2.0f, 0.0f, 0.0f);
        b_y = (d == 2) ? vec3<T>(0.0f, 2.0f, 0.0f) : vec3<T>(0.0f, 0.0f, 2.0f);
    }

    slice(std::array<float, 9> orientation) {
        b_x.x = orientation[0];
        b_x.y = orientation[1];
        b_x.z = orientation[2];

        b_y.x = orientation[3];
        b_y.y = orientation[4];
        b_y.z = orientation[5];

        base.x = orientation[6];
        base.y = orientation[7];
        base.z = orientation[8];
    }

    vec<3_D, T> base; // base point of the sice
    vec<3_D, T> b_x;  // 'basis vector' first direction
    vec<3_D, T> b_y;  // 'basis vector' second direction
};

/** Rotate a 3D vector */
template <typename T>
auto rotate(vec<3_D, T> v, vec<3_D, T> normal, T angle) {
    return glm::rotate(v, angle, normal);
}

/** Rotate a 2D vector */
template <typename T>
auto rotate(vec<2_D, T> v, [[maybe_unused]] vec<2_D, T> normal, T angle) {
    (void)normal;
    return glm::rotate(v, angle);
}

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

/** Check whether an integral value is a power of two. */
template <typename T,
          typename = typename std::enable_if<std::is_integral<T>::value>::type>
bool is_power_of_two(T x) {
    // See Hacker's delight page 11
    return (x & (x - 1)) == 0;
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

/** Checks whether a *voxel-coordinate* is a valid index (is 'inside') the
 * volume. */
template <dimension D, typename T>
bool valid_index(vec<D, int> a, volume<D, T> vol) {
    for (int dim = 0; dim < D; ++dim) {
        if (a[dim] < 0 || a[dim] >= (T)vol.voxels()[dim]) {
            return false;
        }
    }
    return true;
}

/** Checks whether a *voxel-coordinate* vector lies in the box defined by a
 * volume. This is
 * open by one side, the upper side in each axis. */
template <dimension D, typename T>
bool inside(vec<D, T> a, volume<D, T> vol) {
    for (int dim = 0; dim < D; ++dim) {
        if (a[dim] < epsilon<T> ||
            a[dim] > (T)vol.voxels()[dim] - epsilon<T>) {
            return false;
        }
    }
    return true;
}

/** Checks whether a *voxel-coordinate vector* lies in the box defined by a
 * volume. This is
 * open by one side, the upper side in each axis. */
template <dimension D, typename T>
bool inside_margin(vec<D, T> a, volume<D, T> vol, T margin) {
    for (int dim = 0; dim < D; ++dim) {
        if (a[dim] < -(margin + epsilon<T>) ||
            a[dim] >= (T)vol.voxels()[dim] + margin + epsilon<T>)
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
            if (a[d] < -epsilon<T> || a[d] >= sides[d] + epsilon<T>) {
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
vec<D, T> to_voxel(vec<D, T> x, volume<D, T> v) {
    x -= v.origin();
    x /= v.physical_lengths();
    x *= v.voxels();
    return x;
}

template <dimension D, typename T>
optional<line<D, T>> truncate_to_volume(ray<D, T> ray, volume<D, T> v) {
    // need line plane intersection, because the box is axis aligned (AABB) we
    // need to do a ray/AABB intersection.
    auto origin = aabb_intersection<D, T>(
        ray.source, ray.detector, v.physical_lengths(), vec<D, T>(v.origin()));
    if (!origin) {
        return optional<line<D, T>>();
    }

    auto voxel_source = to_voxel<D, T>(ray.source, v);
    auto voxel_detector = to_voxel<D, T>(ray.detector, v);
    auto voxel_origin = to_voxel<D, T>(origin.value().first, v);

    // have to do this in terms of voxel delta
    auto steps = (T)((int)distance(voxel_source, voxel_origin) + 1);
    auto delta = normalize(voxel_detector - voxel_source);
    auto start = voxel_source + steps * delta;

    return optional<line<D, T>>{line<D, T>{start, delta}};
}

/**
 * for a line and bounds, return a pair of vectors consisting
 * of the two intersection points with the bounds
 */
template <dimension D, typename T>
optional<std::pair<vec<D, T>, vec<D, T>>>
intersect_bounds(line<D, T> l, std::array<vec2<int>, D> bounds) {
    vec<D, T> bounds_origin{bounds[0][0], bounds[1][0], bounds[2][0]};

    vec<D, T> bounds_sides =
        vec<D, T>{bounds[0][1] + 1, bounds[1][1] + 1, bounds[2][1] + 1} - bounds_origin;

    auto center = bounds_origin + (T)0.5 * bounds_sides;
    auto for_sure_far_enough = (T)2.0 * distance(l.origin, center); 
    auto b = l.origin + for_sure_far_enough * l.delta;
    auto result = aabb_intersection<D, T>(
        l.origin, b,
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
void interpolate(vec<D, T> a, volume<D, T> v,
                 std::vector<matrix_element<T>>& queue) {
    // First we see what cell corner we are closest to
    vec<D, int> b = round(a);

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

        bool flag = false;
        // if any of the coordinates are outside, we continue
        for (int d = 0; d < D; ++d) {
            if (cell[d] < 0 || cell[d] >= v.voxels()[d]) {
                flag = true;
                break;
            }
        }
        if (flag) {
            continue;
        }

        int index = v.index_by_vector(cell);
        auto value = product<D, T>(vec<D, T>((T)1) - abs(a - cell_center));
        queue.push_back({index, value});
    }
}

template <tomo::dimension D, typename T>
vec<D - 1, T> restrict(vec<D, T> x, int skip) {
    vec<D - 1, T> reduced_point;
    int k = 0;
    for (int d = 0; d < D; ++d) {
        if (d == skip)
            continue;
        reduced_point[k++] = x[d];
    }
    return reduced_point;
}

template <tomo::dimension D, typename T>
vec<D, T> extend(vec<D - 1, T> x, int axis, T value) {
    vec<D, T> point;
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

/** Return the center of the volume */
template <tomo::dimension D, typename T>
vec<D, T> volume_center(tomo::volume<D, T> v) {
    return v.origin() + (T)0.5 * v.physical_lengths();
}

/** Return the optional intersection point of a ray with a plane */
template <typename T>
optional<vec3<T>> ray_plane_intersection(vec3<T> source, vec3<T> direction,
                                         vec3<T> base, vec3<T> normal) {
    auto alpha = dot(normal, direction);
    if (abs(alpha) < 0.001f) {
        // error: ray and plane are parallel
        return {};
    }
    auto offset = dot((base - source), normal) / alpha;
    return {source + offset * direction};
}

} // namespace math
} // namespace tomo
