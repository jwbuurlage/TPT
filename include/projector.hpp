#pragma once

#include <type_traits>
#include <vector>

#include "common.hpp"
#include "geometry.hpp"
#include "logging.hpp"
#include "math.hpp"

namespace tomo {
namespace dim {

/**
 * A discrete integration method (DIM) for the ray inside a volume. Also called
 * a 'projector'
 */
template <dimension D, typename T, typename Iterator>
class base {
  public:
    using problem_dimension = std::integral_constant<dimension, D>;

    /** Construct a projector for a scanning volume. */
    base(volume<D> vol) : volume_(vol) {}

    /** Obtain an iterator to the first voxel that is 'hit'. */
    Iterator begin() { return begin_(); }

    /** Obtain an iterator beyond the final voxel that is 'hit'. */
    Iterator end() { return end_(); }

    /**
     * A container that wraps the `begin` and `end` functions of a projector.
     *
     * This allows for code such as:
     * ```{.cpp}
     * for (auto l : g)
     *      for (auto element : proj(l)) {
     *        // do something with the element
     *      }
     * ```
     */
    class element_container {
      public:
        element_container(base& p) : p_(p) {}
        inline Iterator begin() { return p_.begin(); }
        inline Iterator end() { return p_.end(); }

      private:
        base& p_;
    };

    /** Returns an iterable container for line `line` with this DIM. */
    element_container operator()(math::ray<D, T> incoming_ray) {
        // truncate to volume here first, then reset
        auto truncated_line = truncate_to_volume(incoming_ray, volume_);
        this->clear_();
        if (truncated_line) {
            this->reset_(truncated_line.value());
        }
        return element_container(*this);
    }

    /** Obtain the current line of the DIM. */
    math::line<D, T> get_line() const { return line_; }

    /** Obtain a reference to the scanned volume. */
    volume<D> get_volume() const { return volume_; }

  protected:
    virtual ~base() = default;

    volume<D> volume_;
    math::line<D, T> line_;

  private:
    /* This idiom (keeping virtual functions private) is called NVI */
    virtual Iterator begin_() = 0;
    virtual Iterator end_() = 0;
    virtual void reset_(math::line<D, T> line) = 0;
    virtual void clear_() = 0;

    //TODO: move vector<matrix_elem> queue up here, so we dont need CRT anymore
};

} // namespace dim
} // namespace tomo
