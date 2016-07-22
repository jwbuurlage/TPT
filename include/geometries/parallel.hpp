#pragma once

#include <limits>

#include "../common.hpp"
#include "../geometry.hpp"
#include "../math.hpp"
#include "../volume.hpp"

namespace tomo {

template <dimension Dimension, typename T>
class parallel_line_iterator;

template <dimension Dimension, typename T = default_scalar_type>
class parallel_geometry
    : public geometry<Dimension, T, parallel_line_iterator<Dimension, T>> {};

template <typename T>
class parallel_geometry<2_D, T>
    : public geometry<2_D, T, parallel_line_iterator<2_D, T>> {
  public:
    using iterator_type = parallel_line_iterator<2_D, T>;

    parallel_geometry(int angle_count, int detector_count,
                      const volume<2_D>& volume)
        : geometry<2_D, T, iterator_type>(angle_count * detector_count),
          volume_(volume) {
        auto angle_step = math::pi / angle_count;
        for (T angle = 0.0; angle < math::pi; angle += angle_step) {
            angles_.push_back(angle);
        }

        auto detector_step = (T)volume.y() / detector_count;
        for (int detector = 0; detector < detector_count; detector++) {
            detectors_.push_back((detector - (detector_count - 1) * 0.5) *
                                 detector_step);
        }

        this->dimensions_ = {detector_count, angle_count};
    }

    iterator_type begin() const override { return iterator_type(0, *this); }
    iterator_type end() const override {
        return iterator_type(this->line_count_, *this);
    }

    size_t detector_count() const { return detectors_.size(); }
    size_t angle_count() const { return angles_.size(); }

    const std::vector<T>& angles() const { return angles_; }
    const std::vector<T>& detectors() const { return detectors_; }

    const volume<2_D>& get_volume() const { return volume_; }

  private:
    std::vector<T> angles_;
    std::vector<T> detectors_;
    volume<2_D> volume_;
};

template <dimension Dimension, typename T>
class parallel_line_iterator
    : public std::iterator<std::forward_iterator_tag, line<Dimension, T>> {};

template <typename T>
class parallel_line_iterator<2_D, T>
    : public std::iterator<std::forward_iterator_tag, line<2_D, T>> {
  public:
    parallel_line_iterator(int i, const parallel_geometry<2_D, T>& geometry)
        : i_(i), geometry_(geometry) {}

    parallel_line_iterator(const parallel_line_iterator& other)
        : i_(other.i_), geometry_(other.geometry_) {}

    parallel_line_iterator& operator++(int) {
        auto current = *this;
        ++(*this);
        return current;
    }

    bool operator==(const parallel_line_iterator& other) const {
        return i_ == other.i_;
    }

    bool operator!=(const parallel_line_iterator& other) const {
        return !(*this == other);
    }

    line<2_D, T> operator*() {
        auto current_detector =
            geometry_.detectors()[i_ % geometry_.detector_count()];
        auto current_angle =
            geometry_.angles()[i_ / geometry_.detector_count()];

        // some performance can be gained here by *not* shifting with image
        // center, and maybe we even want to cache these results somehow
        auto vol = geometry_.get_volume();
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

        auto end = origin + 3.0 * vol.x() * delta;

        T min_distance = std::numeric_limits<T>::max();
        math::vec2<T> best_point;
        for (auto& line_segment : lines) {
            auto intersection_point = math::intersection<T>(
                origin, end, line_segment[0], line_segment[1]);
            if (intersection_point != math::vec2<T>(0, 0)) {
                auto dist =
                    math::distance<2_D, double>(origin, intersection_point);
                if (dist < min_distance) {
                    best_point = intersection_point;
                    min_distance = dist;
                }
            }
        }

        return {best_point, delta};
    }

    parallel_line_iterator& operator++() {
        i_++;
        return *this;
    }

  private:
    int i_;
    const parallel_geometry<2_D, T>& geometry_;
};

} // namespace tomo
