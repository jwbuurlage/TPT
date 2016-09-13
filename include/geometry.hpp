#pragma once

#include <array>
#include <cassert>
#include <vector>

#include "math.hpp"

namespace tomo {

template <dimension D, typename T = default_scalar_type>
struct line {
    line() = default;

    line(math::vec<D, T> origin_in, math::vec<D, T> delta_in) : origin(origin_in), delta(delta_in) {}

    line(std::initializer_list<math::vec<D, T>> vectors)
        : line(*vectors.begin(), *(vectors.begin() + 1)) {
        // "trying to construct line with wrong number of vectors"
        assert(vectors.size() == 2);
    }

    // the direction vector
    math::vec<D, T> origin;

    // origin vector, first intersection of line with volume
    math::vec<D, T> delta;
};

// => Q: how do we deal with accelerators / distributed in this context
// Possible A: we still generate the begin/end points of the line on the
// CPU,
// but the 'kernel iterators' are implemented on the accelerator
template <dimension D, typename T, typename Derived>
class geometry {
  public:
    class geometry_iterator
        : public std::iterator<std::forward_iterator_tag, line<D, T>> {
      public:
        geometry_iterator(int i, const Derived& geometry)
            : i_(i), geometry_(geometry) {}

        geometry_iterator(const geometry_iterator& other)
            : i_(other.i_), geometry_(other.geometry_) {}

        geometry_iterator& operator++(int) {
            auto current = *this;
            ++(*this);
            return current;
        }

        bool operator==(const geometry_iterator& other) const {
            return i_ == other.i_;
        }

        bool operator!=(const geometry_iterator& other) const {
            return !(*this == other);
        }

        line<D, T> operator*() const { return geometry_.get_line(i_); }

        geometry_iterator& operator++() {
            i_++;
            return *this;
        }

      private:
        int i_;
        const Derived& geometry_;
    };

    geometry(int line_count) : line_count_(line_count) {}

    geometry_iterator begin() const {
        return geometry_iterator(0, *((Derived*)this));
    }
    geometry_iterator end() const {
        return geometry_iterator(line_count_, *((Derived*)this));
    }

    int lines() const { return line_count_; }

    math::vec2<int> groups() const { return dimensions_; }

  protected:
    int line_count_;
    math::vec2<int> dimensions_;
};

} // namespace tomo
