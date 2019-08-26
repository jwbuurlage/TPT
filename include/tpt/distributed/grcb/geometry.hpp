#pragma once

#include <array>
#include <utility>

#include "../../geometry.hpp"
#include "../../math.hpp"

#include <boost/geometry.hpp>
namespace bg = boost::geometry;

namespace tpt::grcb {

// Detector coordinates: (0, 0) is center of detector,
// otherwise follows same distances as in 'world' coordinates

/** Get detector polygon */
template <typename T>
bg::model::polygon<math::vec2<T>> detector(geometry::projection<3_D, T> pi) {
    auto hds = (T)0.5 * pi.detector_size;
    return bg::model::polygon<math::vec2<T>>{
        {-hds, math::vec2<T>{-hds[0], hds[1]}, hds, {hds[0], -hds[1]}, -hds},
        {}};
}

/** Project point onto detector. */
template <typename T>
auto project(geometry::projection<3_D, T> pi, vec3<T> x) {
    return math::ray_plane_intersection<T>(
        pi.source_location, x - pi.source_location, pi.detector_location,
        math::cross<T>(pi.detector_tilt[0], pi.detector_tilt[1]));
}

/** Convert from world coordinates to detector coordinates. */
template <typename T>
math::vec2<T> world_to_detector(geometry::projection<3_D, T> pi, vec3<T> x) {
    auto z = x - pi.detector_location;
    auto u = math::dot<3_D, T>(pi.detector_tilt[0], z);
    auto v = math::dot<3_D, T>(pi.detector_tilt[1], z);
    return {u, v};
}

/** Obtain the shadow for a part */
template <typename T>
std::pair<std::array<math::vec2<T>, 8>, bg::model::polygon<math::vec2<T>>>
shadow(geometry::projection<3_D, T> pi, cube<T> v) {
    auto p = std::array<math::vec2<T>, 8>{};
    std::transform(v.begin(), v.end(), p.begin(), [&](auto x) {
        return world_to_detector(pi, project(pi, x).value());
    });
    auto hull = bg::model::polygon<math::vec2<T>>{};
    bg::convex_hull(p, hull);
    return {p, hull};
}

/** Compute the overlap of the shadows of two parts for a projection. */
template <typename T>
T overlap(geometry::projection<3_D, T> pi, cube<T> vL, cube<T> vR) {
    auto [pL, shadowL] = shadow(pi, vL);
    auto [pR, shadowR] = shadow(pi, vR);

    auto det = detector(pi);

    auto overlap = std::vector<bg::model::polygon<math::vec2<T>>>{};
    bg::intersection(shadowL, shadowR, overlap);

    /** Even though there should be a single intersection since the hulls (and
     detector are convex, we need this loop since there seem to be some line
     segments in the output */
    auto ans = 0.0;
    for (auto& ovl : overlap) {
        auto restricted_overlap =
            std::vector<bg::model::polygon<math::vec2<T>>>{};
        bg::intersection(ovl, det, restricted_overlap);
        for (auto& rovl : restricted_overlap) {
            ans += bg::area(rovl);
        }
    }

    return ans;
}

} // namespace tpt::grcb
