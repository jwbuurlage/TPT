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
    projector(line<D, T> line, volume<D> vol) : line_(line), volume_(vol) {}

    Iterator begin() { return Iterator(*this); }
    Iterator end() { return Iterator(*this, true); }

    line<D, T> get_line() const { return line_; }
    volume<D> get_volume() const { return volume_; }

  private:
    line<D, T> line_;
    volume<D> volume_;
};

template <dimension D, typename T = default_scalar_type>
class linear_projector_iterator
    : public std::iterator<std::forward_iterator_tag, line<D, T>> {
  public:
    using Base = projector<D, T, linear_projector_iterator<D, T>>;
    using vecD = typename math::vec<D, T>::type;

    linear_projector_iterator(const Base& projector)
        : projector_(projector) {
        current_point_ = projector_.get_line().origin;

        // i think this is specific to 2d
        queue_.reserve((int)(math::sqrt(2.0) * 4 * projector_.get_volume()[0]));

        fill_queue();
    }

    // used for end
    linear_projector_iterator(const Base& projector, bool)
        : projector_(projector) {}

    linear_projector_iterator(const linear_projector_iterator& other)
        : projector_(other.projector_), current_point_(other.current_point_),
          queue_(other.queue_) {}

    linear_projector_iterator& operator++(int) {
        auto current = *this;
        ++(*this);
        return current;
    }

    bool operator==(const linear_projector_iterator& other) const {
        return (i_ == queue_.size()) == (other.i_ == other.queue_.size());
    }

    bool operator!=(const linear_projector_iterator& other) const {
        return !(*this == other);
    }

    math::matrix_element<T> operator*() { return queue_[i_]; }

    linear_projector_iterator& operator++() {
        i_++;
        return *this;
    }

    void fill_queue() {
        if (!math::inside<D, T>(current_point_, projector_.get_volume())) {
            TOMO_LOG_VAR(current_point_.x);
            TOMO_LOG_VAR(current_point_.y);
            TOMO_LOG_VAR("queue");
        }
        while (math::inside<D, T>(current_point_, projector_.get_volume())) {
            auto new_points = math::interpolate<D, T>(current_point_,
                                                      projector_.get_volume());
            queue_.insert(queue_.end(), new_points.begin(), new_points.end());
            current_point_ += projector_.get_line().delta;
        }
    }

  private:
    const Base& projector_;
    size_t i_ = 0;

    vecD current_point_;

    std::vector<math::matrix_element<T>> queue_;
};

template <dimension D, typename T = default_scalar_type>
using linear_projector = projector<D, T, linear_projector_iterator<D, T>>;

} // namespace tomo
