#pragma once

#include <limits>

#include "../geometries/trajectory.hpp"
#include "../math.hpp"

namespace tomo {
namespace distributed {

struct shadow {
    math::vec2<int> min_pt;
    math::vec2<int> max_pt;
};

inline auto empty(shadow s) {
    return (s.max_pt.x < s.min_pt.x) || (s.max_pt.y < s.min_pt.y);
}

inline auto pixels(shadow s) {
    return (s.max_pt.x - s.min_pt.x + 1) * (s.max_pt.y - s.min_pt.y + 1);
}

template <typename T>
math::vec<2_D, int> compute_pixel_intersection(math::vec3<T> x,
                                               geometry::projection<3_D, T> p) {
    auto geometric_point = math::ray_plane_intersection<T>(
        p.source_location, x - p.source_location, p.detector_location,
        math::cross<T>(p.detector_tilt[0], p.detector_tilt[1]));
    assert(geometric_point);
    auto pt = geometric_point.value();
    // detector location is the *center* of the detector, so we adjust
    auto offset = pt - (p.detector_location -
                        (T)0.5 * (p.detector_size[0] * p.detector_tilt[0] +
                                  p.detector_size[1] * p.detector_tilt[1]));
    auto dx = math::dot<3_D, T>(offset, p.detector_tilt[0]);
    auto dy = math::dot<3_D, T>(offset, p.detector_tilt[1]);
    return {(int)(dx / p.detector_size[0] * p.detector_shape[0] + (T)0.5),
            (int)(dy / p.detector_size[1] * p.detector_shape[1] + (T)0.5)};
}

template <typename T>
shadow compute_shadow(std::vector<math::vec3<T>> xs,
                      geometry::projection<3_D, T> p) {
    shadow s = {
        {std::numeric_limits<int>::max(), std::numeric_limits<int>::max()},
        {-1, -1}};
    for (auto x : xs) {
        auto intersection = compute_pixel_intersection(x, p);
        s.min_pt = math::min(s.min_pt, intersection);
        s.max_pt = math::max(s.max_pt, intersection);
    }
    s.min_pt = math::max(s.min_pt - math::vec2<int>{1, 1}, {0, 0});
    s.max_pt = math::min(s.max_pt + math::vec2<int>{1, 1},
                         p.detector_shape - math::vec2<int>{1, 1});
    return s;
}

template <typename T>
class restricted_geometry : public geometry::base<3, T> {
  public:
    restricted_geometry(geometry::trajectory<3_D, T>& geometry,
                        volume<3_D, T> local_volume)
        : geometry::base<3, T>(geometry.projection_count(), false),
          geometry_(geometry), local_volume_(local_volume),
          pixels_(geometry.projection_count() + 1) {
        project();
        std::transform(shadows_.begin(), shadows_.end(), pixels_.begin() + 1,
                       [&](auto shadow) { return pixels(shadow); });
        std::partial_sum(pixels_.begin(), pixels_.end(), pixels_.begin());
        this->compute_lines_();
    }

    void project() {
        for (int i = 0; i < geometry_.projection_count(); ++i) {
            shadows_.push_back(compute_shadow(local_volume_.corners(),
                                              geometry_.get_projection(i)));
        }
    }

    math::vec<2_D, int> projection_shape(int i) const override {
        auto dx = shadows_[i].max_pt.x + 1 - shadows_[i].min_pt.x;
        auto dy = shadows_[i].max_pt.y + 1 - shadows_[i].min_pt.y;
        if (dx <= 0 || dy <= 0) {
            return {0, 0};
        }
        return {dx, dy};
    }

    math::vec<3_D, T> detector_corner(int i) const override {
        return geometry_.detector_corner(i) +
               (T)shadows_[i].min_pt[0] * geometry_.projection_delta(i)[0] +
               (T)shadows_[i].min_pt[1] * geometry_.projection_delta(i)[1];
    }

    math::vec<3_D, T> source_location(int i) const override {
        return geometry_.source_location(i);
    }

    std::array<math::vec<3_D, T>, 2> projection_delta(int i) const override {
        return geometry_.projection_delta(i);
    }

    auto local_shadow(int i) const { return shadows_[i]; }
    const auto& global_geometry() const { return geometry_; }

  private:
    geometry::trajectory<3_D, T>& geometry_;
    volume<3_D, T> local_volume_;
    std::vector<shadow> shadows_;
    std::vector<int> pixels_;
};

} // namespace distributed
} // namespace tomo
