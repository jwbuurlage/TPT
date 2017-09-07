#pragma once

#include <array>
#include <cassert>
#include <cstdint>
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
    using problem_dimension = std::integral_constant<dimension, D>;
    using value_type = T;

    /** An iterator to a line of the geometry. */
    class pixel_iterator
        : public std::iterator<std::forward_iterator_tag, math::ray<D, T>> {
      public:
        pixel_iterator() {}

        /** Construct the iterator with a projection index and a geometry. */
        pixel_iterator(math::vec<D - 1, int> projection_shape,
                       math::vec<D, T> detector_corner,
                       math::vec<D, T> source_location,
                       std::array<math::vec<D, T>, D - 1> projection_delta,
                       bool parallel = false, math::vec<D - 1, int> pixel = {})
            : projection_shape_(projection_shape),
              detector_corner_(detector_corner),
              source_location_(source_location),
              projection_delta_(projection_delta),
              current_location_(
                  detector_corner_ +
                  (T)0.5 * (projection_delta[0] + projection_delta[1])),
              parallel_(parallel), pixel_(pixel) {}

        /** Copy-construct an iterator. */
        pixel_iterator(const pixel_iterator& other)
            : projection_shape_(other.projection_shape_),
              detector_corner_(other.detector_corner_),
              source_location_(other.source_location_),
              projection_delta_(other.projection_delta_),
              current_location_(other.current_location_),
              parallel_(other.parallel_), pixel_(other.pixel_) {}

        /** Increase the iterator. */
        pixel_iterator& operator++(int) {
            auto current = *this;
            ++(*this);
            return current;
        }

        /** Compare two iterators for inequality. */
        bool operator!=(const pixel_iterator& other) const {
            return !(*this == other);
        }

        /** Compare two iterators for equality. */
        bool operator==(const pixel_iterator& other) const {
            return pixel_ == other.pixel_;
        }

        /** Dereference the iterator to obtain the represented line. */
        math::ray<D, T> operator*() const {
            if (!parallel_) {
                return {source_location_, current_location_};
            } else {
                return {source_location_ + current_location_ - detector_corner_,
                        current_location_};
            }
        }

        /** Increase the iterator. */
        pixel_iterator& operator++() {
            // FIXME check if this gets unrolled
            for (int d = 0; d < D - 2; ++d) {
                pixel_[d] += 1;
                if (pixel_[d] < projection_shape_[d]) {
                    current_location_ += projection_delta_[d];
                    return *this;
                }
                pixel_[d] = 0;
                // reset d-th dimension
                current_location_ -=
                    (T)(projection_shape_[d] - 1) * projection_delta_[d];
            }
            pixel_[D - 2]++;
            current_location_ += projection_delta_[D - 2];
            return *this;
        }

      private:
        math::vec<D - 1, int> projection_shape_ = {};
        math::vec<D, T> detector_corner_ = {};
        math::vec<D, T> source_location_ = {};
        std::array<math::vec<D, T>, D - 1> projection_delta_ = {};
        math::vec<D, T> current_location_ = {};
        bool parallel_ = false;
        math::vec<D - 1, int> pixel_ = {};
    };

    /** An iterator for the projections of the geometry. */
    class projection_iterator
        : public std::iterator<std::forward_iterator_tag, math::ray<D, T>> {
      public:
        /** Construct the iterator with a line index and a geometry. */
        projection_iterator(int proj, const base& geometry)
            : proj_(proj), geometry_(geometry), pixels_(first_()),
              pixels_last_(last_()) {
            this->update_();
        }

        /** Copy-construct an iterator. */
        projection_iterator(const projection_iterator& other)
            : proj_(other.proj_), geometry_(other.geometry_) {}

        /** Increase the iterator. */
        projection_iterator& operator++(int) {
            auto current = *this;
            ++(*this);
            return current;
        }

        /** Compare two iterators for inequality. */
        bool operator!=(const projection_iterator& other) const {
            return !(*this == other);
        }

        /** Compare two iterators for equality. */
        bool operator==(const projection_iterator& other) const {
            return proj_ == other.proj_ && pixels_ == other.pixels_;
        }

        /** Dereference the iterator to obtain the represented line. */
        math::ray<D, T> operator*() const { return *pixels_; }

        /** Increase the iterator. */
        projection_iterator& operator++() {
            ++pixels_;
            this->update_();
            return *this;
        }

      private:
        void update_() {
            while (pixels_ == pixels_last_ &&
                   proj_ < geometry_.projection_count()) {
                proj_++;
                pixels_ = first_();
                pixels_last_ = last_();
            }
        }

        inline pixel_iterator first_() {
            if (proj_ >= geometry_.projection_count()) {
                return pixel_iterator();
            }
            auto shape = geometry_.projection_shape(proj_);
            auto corner = geometry_.detector_corner(proj_);
            auto location = geometry_.source_location(proj_);
            auto delta = geometry_.projection_delta(proj_);
            auto parallel = geometry_.parallel();
            return pixel_iterator(shape, corner, location, delta, parallel);
        }

        inline pixel_iterator last_() {
            if (proj_ >= geometry_.projection_count()) {
                return pixel_iterator();
            }
            auto shape = geometry_.projection_shape(proj_);
            auto corner = geometry_.detector_corner(proj_);
            auto location = geometry_.source_location(proj_);
            auto delta = geometry_.projection_delta(proj_);
            auto parallel = geometry_.parallel();
            math::vec<D - 1, int> outside = {};
            outside[D - 2] = shape[D - 2];
            return pixel_iterator(shape, corner, location, delta, parallel,
                                  outside);
        }

        int proj_ = 0;
        const base& geometry_;
        pixel_iterator pixels_;
        pixel_iterator pixels_last_;
    };

    /** Construct the geometry with a given number of lines. */
    base(int projection_count, bool parallel = false)
        : projection_count_(projection_count), parallel_(parallel) {}

    virtual ~base() = default;

    /** Obtain an iterator to the first element of the geometry. */
    projection_iterator begin() const { return projection_iterator(0, *this); }

    /** Obtain an iterator beyond the last element of the geometry. */
    projection_iterator end() const {
        return projection_iterator(projection_count_, *this);
    }

    /** Obtain the number of lines */
    auto lines() const { return line_count_; }

    /** Obtain the projection count */
    int projection_count() const { return projection_count_; }

    bool parallel() const { return parallel_; }

    virtual math::vec<D - 1, int> projection_shape(int i) const = 0;
    virtual math::vec<D, T> detector_corner(int i) const = 0;
    virtual math::vec<D, T> source_location(int i) const = 0;

    virtual std::array<math::vec<D, T>, D - 1>
    projection_delta(int i) const = 0;
    int offset(int idx) const { return offsets_[idx]; }

  protected:
    std::vector<int> offsets_;

    void compute_offsets_() {
        offsets_.resize(projection_count());
        offsets_[0] = 0;
        std::iota(offsets_.begin() + 1, offsets_.end(), 0);
        std::transform(
            offsets_.begin() + 1, offsets_.end(), offsets_.begin() + 1,
            [&](int i) { return math::reduce<D - 1>(projection_shape(i)); });
        std::partial_sum(offsets_.begin() + 1, offsets_.end(),
                         offsets_.begin() + 1);
    }

    void compute_lines_() {
        for (auto i = 0; i < projection_count_; ++i) {
            line_count_ += math::reduce<D - 1>(this->projection_shape(i));
        }
        this->compute_offsets_();
    }

    int projection_count_ = 0;
    uint64_t line_count_ = 0;
    bool parallel_ = false;
};

} // namespace geometry
} // namespace tomo
