#pragma once

#include <limits>

#include "../common.hpp"
#include "../geometry.hpp"
#include "../math.hpp"
#include "../volume.hpp"

namespace tomo {

template <typename T>
T detector_location(int detector, int detector_count, T detector_step,
                    const volume<2_D>&) {
    return (detector - (detector_count - 1) * 0.5) * detector_step;
}

template <typename T>
typename math::vec<2_D, T>::type
detector_location(int detector, int detector_count, T detector_step,
                  const volume<3_D>&) {
    auto detector_x = detector % detector_count;
    auto detector_y = detector / detector_count;

    return {(detector_x - (detector_count - 1) * 0.5) * detector_step,
            (detector_y - (detector_count - 1) * 0.5) * detector_step};
}

template <dimension D, typename T>
class parallel_geometry : public geometry<D, T, parallel_geometry<D, T>> {
  public:
    using position = math::vec<D - 1, T>;

    parallel_geometry(int angle_count, int detector_count,
                      const volume<D>& volume)
        : geometry<D, T, parallel_geometry<D, T>>(
              angle_count * math::pow(detector_count, D - 1)),
          volume_(volume) {
        auto angle_step = math::pi<T> / angle_count;
        for (T angle = 0.0; angle < math::pi<T>; angle += angle_step) {
            angles_.push_back(angle);
        }

        int total_detector_count = math::pow(detector_count, D - 1);
        // FIXME this is only for equilateral volume
        auto detector_step = volume_.y() / (T)detector_count;
        for (int detector = 0; detector < total_detector_count; detector++) {
            detectors_.push_back(detector_location<T>(detector, detector_count,
                                                      detector_step, volume_));
        }

        this->dimensions_ = {detector_count, angle_count};
    }

    size_t detector_count() const { return detectors_.size(); }
    size_t angle_count() const { return angles_.size(); }

    const std::vector<T>& angles() const { return angles_; }
    const std::vector<position>& detectors() const { return detectors_; }

    const volume<D>& get_volume() const { return volume_; }

    line<D, T> get_line(int i) const {
        return compute_line(detectors_[i % detector_count()],
                            angles_[i / detector_count()], volume_);
    }

  private:
    std::vector<T> angles_;
    std::vector<position> detectors_;
    volume<D> volume_;
};

// template <dimension D, typename T>
// inline line<D, T> compute_line(T, T, const volume<D>&) {
//    static_assert(false,
//                  "computing lines for `D /= 2` dimensions not supported");
//}

template <typename T>
inline line<2_D, T> compute_line(T current_detector, T current_angle,
                                 const volume<2_D>& vol) {
    // some performance can be gained here by *not* shifting with image
    // center, and maybe we even want to cache these results somehow
    auto origin = math::vec2<T>(-vol.x(), current_detector);

    auto c = math::cos(-current_angle);
    auto s = math::sin(-current_angle);

    origin = math::vec2<T>(c * origin[0] - s * origin[1],
                           s * origin[0] + c * origin[1]);

    auto image_center = math::vec2<T>(0.5 * vol.x(), 0.5 * vol.y());
    origin += image_center;

    auto delta = math::vec2<T>(c, s);

    // we now need to intersect the line current, that has an angle
    // current_angle_ with our box
    // we can intersect with all 'three top lines' (a, b, c) and choose
    // the intersection point with smallest distance

    //  _->
    // /
    //   ----b----
    //   |       |
    // + |		 | -
    // + a		 c -
    // + | 		 | -
    //   |       |
    //	 ----d---- /
    //            /
    //         <--

    static const std::array<std::array<math::vec2<T>, 2>, 3> lines = {
        std::array<math::vec2<T>, 2>{math::vec2<T>(0, 0),
                                     math::vec2<T>(0, vol.y())},
        std::array<math::vec2<T>, 2>{math::vec2<T>(0, vol.y()),
                                     math::vec2<T>(vol.x(), vol.y())},
        std::array<math::vec2<T>, 2>{math::vec2<T>(vol.x(), vol.y()),
                                     math::vec2<T>(vol.x(), 0)}};

    auto end = origin + (T)3.0 * vol.x() * delta;

    T min_distance = std::numeric_limits<T>::max();
    math::vec2<T> best_point;
    for (auto& line_segment : lines) {
        auto intersection_point = math::intersection<T>(
            origin, end, line_segment[0], line_segment[1]);
        if (intersection_point != math::vec2<T>(0, 0)) {
            auto dist = math::distance<2_D, T>(origin, intersection_point);
            if (dist < min_distance) {
                best_point = intersection_point;
                min_distance = dist;
            }
        }
    }

    return {best_point + (T)0.5 * delta, delta};
}

template <typename T>
inline line<3_D, T> compute_line(math::vec2<T> current_detector,
                                 T current_angle, const volume<3_D>& vol) {
    // strategy: only consider current detector x, and ignore y, only add it at
    // the end
    auto volume_slice = volume<2_D>(vol.x(), vol.y());
    auto line_2d =
        compute_line(current_detector.x, current_angle, volume_slice);

    auto origin = math::vec3<T>(line_2d.origin.x, line_2d.origin.y,
                                current_detector.y + 0.5 * vol.z());
    auto delta = math::vec3<T>(line_2d.delta.x, line_2d.delta.y, 0.0);

    return {origin, delta};
}

} // namespace tomo
