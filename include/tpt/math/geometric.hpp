#pragma once

#include "vector.hpp"
#include "vector_operations.hpp"

#include "../geometry.hpp"

namespace tpt {
using namespace tpt::literals;

namespace math {

template <typename T>
vec<1_D, int> project_hyxel(vec2<T> x, geometry::projection<2_D, T> p) {

    /*
     * a --------*-------------- b
     *            \
     *             \
     *              x
     *               \
     *             source
     */

    auto offset = (T)0.5 * p.detector_size[0] * p.detector_tilt[0];
    auto d = p.detector_location;
    auto a = d - offset;
    auto b = d + offset;
    auto s = p.source_location;

    auto detector_point = line_intersection<T>(a, b, s, x);
    assert(detector_point);

    auto pixel =
        (int)(p.detector_shape[0] *
              math::dot<2_D, T>(detector_point.value() - a, p.detector_tilt[0]) /
              p.detector_size[0]);

    return vec1<int>(pixel);
}

template <typename T>
vec<2_D, int> project_hyxel(vec3<T> x, geometry::projection<3_D, T> p) {
    auto geometric_point = ray_plane_intersection<T>(
        p.source_location, x - p.source_location, p.detector_location,
        math::cross<T>(p.detector_tilt[0], p.detector_tilt[1]));
    assert(geometric_point);
    auto pt = geometric_point.value();
    // detector location is the *center* of the detector, so we adjust
    auto offset = pt - (p.detector_location -
                        (T)0.5 * (p.detector_size[0] * p.detector_tilt[0] +
                                  p.detector_size[1] * p.detector_tilt[1]));
    auto dx = dot<3_D, T>(offset, p.detector_tilt[0]);
    auto dy = dot<3_D, T>(offset, p.detector_tilt[1]);
    return {(int)(dx / p.detector_size[0] * p.detector_shape[0] + (T)0.5),
            (int)(dy / p.detector_size[1] * p.detector_shape[1] + (T)0.5)};
}

template <typename T>
vec<1_D, int> parallel_project_hyxel(vec2<T> x,
                                     geometry::projection<2_D, T> p) {
    std::cerr << "`parallel_project_hyxel` NOT IMPLEMENTED\n";
    (void)x;
    (void)p;
    return vec1<int>(0);
}

template <typename T>
vec<2_D, int> parallel_project_hyxel(vec3<T> x,
                                     geometry::projection<3_D, T> p) {
    std::cerr << "`parallel_project_hyxel` NOT IMPLEMENTED\n";
    (void)x;
    (void)p;
    return vec2<int>(0);
}

} // namespace math
} // namespace tpt
