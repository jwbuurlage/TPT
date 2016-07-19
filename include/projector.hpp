#pragma once

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
        : projector_(projector), finished_(false) {
        current_point_ = projector_.get_line().origin;
        fill_queue();
        next_element();
    }

    // used for end
    linear_projector_iterator(const Base& projector, bool finished)
        : projector_(projector), finished_(finished) {}

    linear_projector_iterator(const linear_projector_iterator& other)
        : projector_(other.projector_), current_point_(other.current_point_),
          finished_(other.finished_), queue_(other.queue_) {}

    linear_projector_iterator& operator++(int) {
        auto current = *this;
        ++(*this);
        return current;
    }

    bool operator==(const linear_projector_iterator& other) const {
        return finished_ == other.finished_;
    }

    bool operator!=(const linear_projector_iterator& other) const {
        return !(*this == other);
    }

    math::matrix_element<T> operator*() {
        if (finished_)
            return math::matrix_element<T>{-1, (T)0};
        return current_element_;
    }

    linear_projector_iterator& operator++() {
        if (queue_.empty()) {
            current_point_ += projector_.get_line().delta;
            if (!math::inside<D, T>(current_point_, projector_.get_volume())) {
                finished_ = true;
                return *this;
            }
            fill_queue();
        }

        next_element();
        return *this;
    }

    void fill_queue() {
        queue_ =
            math::interpolate<D, T>(current_point_, projector_.get_volume());
    }

    void next_element() {
        if (queue_.empty())
            return;

        current_element_ = queue_.back();
        queue_.pop_back();
    }

  private:
    const Base& projector_;

    vecD current_point_;
    math::matrix_element<T> current_element_ =
        math::matrix_element<T>{-1, (T)0};

    bool finished_;
    std::vector<math::matrix_element<T>> queue_;
};

template <dimension D, typename T = default_scalar_type>
using linear_projector = projector<D, T, linear_projector_iterator<D, T>>;

} // namespace tomo
