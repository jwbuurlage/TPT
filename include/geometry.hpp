#pragma once

#include <array>
#include <cassert>
#include <vector>
#include <type_traits>

#include "math.hpp"

namespace tomo {

/**
 * A line inside a volume.
 *
 * \tparam D the dimension of the volume.
 * \tparam T the scalar type to use
 */
template <dimension D, typename T = default_scalar_type>
struct line {
    /** Construct a line using defaults. */
    line() = default;

    /** Construct a line from two vectors. */
    line(math::vec<D, T> origin_in, math::vec<D, T> delta_in)
        : origin(origin_in), delta(delta_in) {}

    /** Construct a line from an initializer list. */
    line(std::initializer_list<math::vec<D, T>> vectors)
        : line(*vectors.begin(), *(vectors.begin() + 1)) {
        // "trying to construct line with wrong number of vectors"
        assert(vectors.size() == 2);
    }

    math::vec<D, T> origin; //> first intersection of line with a volume
    math::vec<D, T> delta;  //> the direction vector
};

namespace geometry {

/**
 * A description of the experimental setup.
 *
 * The geometry is a glorified container, which describes all the rays that are
 * used for the object reconstruction.
 *
 * \tparam D the dimension of the volume.
 * \tparam T the scalar type to use
 * \tparam Derived the type of the actual geometry, used for CRT inside
 * `geomety_iterator`.
 */
template <dimension D, typename T, typename Derived>
class base {
  public:
    using problem_dimension = std::integral_constant<dimension, D>;

    /** An iterator to a line of the geometry. */
    class geometry_iterator
        : public std::iterator<std::forward_iterator_tag, line<D, T>> {
      public:

        /** Construct the iterator with a line index and a geometry. */
        geometry_iterator(int i, const Derived& geometry)
            : i_(i), geometry_(geometry) {}

        /** Copy-construct an iterator. */
        geometry_iterator(const geometry_iterator& other)
            : i_(other.i_), geometry_(other.geometry_) {}

        /** Increase the iterator. */
        geometry_iterator& operator++(int) {
            auto current = *this;
            ++(*this);
            return current;
        }

        /** Compare two iterators for equality. */
        bool operator==(const geometry_iterator& other) const {
            return i_ == other.i_;
        }

        /** Compare two iterators for inequality. */
        bool operator!=(const geometry_iterator& other) const {
            return !(*this == other);
        }

        /** Dereference the iterator to obtain the represented line. */
        line<D, T> operator*() const { return geometry_.get_line(i_); }

        /** Increase the iterator. */
        geometry_iterator& operator++() {
            i_++;
            return *this;
        }

      private:
        int i_;
        const Derived& geometry_;
    };

    /** Construct the geometry with a given number of lines. */
    base(int line_count) : line_count_(line_count) {}

    /** Obtain an iterator to the first element of the geometry. */
    geometry_iterator begin() const {
        return geometry_iterator(0, *((Derived*)this));
    }
    /** Obtain an iterator beyond the last element of the geometry. */
    geometry_iterator end() const {
        return geometry_iterator(line_count_, *((Derived*)this));
    }

    /** Obtain the number of lines */
    int lines() const { return line_count_; }

    /**
     * Obtain the grouping of a geometry.
     *
     * Lines within the geometry can be grouped together. For example, lines may
     * correspond to the same view or have roughly the same angle. This can be
     * used by algorithms, or make for easier visualization of the measurement
     * data.
     */
    math::vec2<int> groups() const { return dimensions_; }

  protected:
    int line_count_;
    math::vec2<int> dimensions_;
};

} // namespace geometry
} // namespace tomo
