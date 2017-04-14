#pragma once

#include <array>
#include <cassert>
#include <type_traits>
#include <vector>

#include "math.hpp"

namespace tomo {
namespace geometry {

/**
 * A description of the geometry of the acquisition.
 *
 * The geometry is a glorified container, which describes all the rays that are
 * used for the object reconstruction.
 *
 * \tparam D the dimension of the volume.
 * \tparam T the scalar type to use
 */
template <dimension D, typename T>
class base {
  public:
    /** An iterator to a line of the geometry. */
    class geometry_iterator
        : public std::iterator<std::forward_iterator_tag, math::ray<D, T>> {
      public:
        /** Construct the iterator with a line index and a geometry. */
        geometry_iterator(int i, const base& geometry)
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
        math::ray<D, T> operator*() const { return geometry_.get_line(i_); }

        /** Increase the iterator. */
        geometry_iterator& operator++() {
            i_++;
            return *this;
        }

      private:
        int i_;
        const base& geometry_;
    };

    using problem_dimension = std::integral_constant<dimension, D>;
    using value_type = T;

    /** Construct the geometry with a given number of lines. */
    base(int projection_count, math::vec<D - 1, int> detector_shape)
        : projection_count_(projection_count), detector_shape_(detector_shape),
          line_count_(projection_count * math::reduce<D - 1>(detector_shape)),
          detector_pixel_count_(math::product<D - 1, int>(detector_shape)) {}

    virtual ~base() = default;

    /** Obtain an iterator to the first element of the geometry. */
    geometry_iterator begin() const { return geometry_iterator(0, *this); }

    /** Obtain an iterator beyond the last element of the geometry. */
    geometry_iterator end() const {
        return geometry_iterator(line_count_, *this);
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
    math::vec<D - 1, int> detector_shape() const { return detector_shape_; }
    void set_detector_shape(math::vec<D - 1, int> detector_shape) {
        detector_shape_ = detector_shape;
        recompute_lines_();
    }

    /* Obtain the i-th line */
    virtual math::ray<D, T> get_line(int i) const = 0;

    int projection_count() const { return projection_count_; }
    int detector_pixel_count() const { return detector_pixel_count_; }

    virtual void set_projections(int projection_count) {
        projection_count_ = projection_count;
        recompute_lines_();
    }

  protected:
    void recompute_lines_() {
        detector_pixel_count_ = math::product<D - 1, int>(detector_shape_);
        this->line_count_ =
            this->projection_count_ * this->detector_pixel_count_;
    }

    int projection_count_;
    math::vec<D - 1, int> detector_shape_;
    int line_count_;
    int detector_pixel_count_;
};

} // namespace geometry
} // namespace tomo
