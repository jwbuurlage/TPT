#pragma once

#include <vector>

#include "common.hpp"
#include "geometry.hpp"
#include "logging.hpp"
#include "math.hpp"

namespace tomo {

template <dimension D, typename T, typename Iterator>
class projector {
  public:
    projector(volume<D> vol) : volume_(vol), line_({{0, 0}, {0, 0}}) {
        queue_.reserve((int)(math::sqrt2 * 4 * volume_[0]));
    }

    Iterator begin() { return queue_.begin(); }
    Iterator end() { return queue_.end(); }

    void reset(line<D, T> line) {
        queue_.clear();
        auto current_point = line.origin;

        while (math::inside<D, T>(current_point, volume_)) {
            math::interpolate<D, T>(current_point, volume_, queue_);
            current_point += line.delta;
        }
    }

    line<D, T> get_line() const { return line_; }
    volume<D> get_volume() const { return volume_; }

  private:
    volume<D> volume_;
    line<D, T> line_;

    std::vector<math::matrix_element<T>> queue_;
};

template <dimension D, typename T = default_scalar_type>
using linear_projector =
    projector<D, T, typename std::vector<math::matrix_element<T>>::iterator>;

} // namespace tomo
