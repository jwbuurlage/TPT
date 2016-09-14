#pragma once

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
    /** Construct a projector for a scanning volume. */
    base(volume<D> vol) : volume_(vol) {}

    /** Obtain an iterator to the first voxel that is 'touched'. */
    Iterator begin() { return begin_(); }

    /** Obtain an iterator beyond the final voxel that is 'touched'. */
    Iterator end() { return end_(); }

    /**
     * A container that wraps the `begin` and `end` functions of a projector.
     *
     * This allows for code such as:
     * ```{.cpp}
     * for (auto l : g)
     *      for (auto touch : proj(l)) {
     *        // do something with the touch
     *      }
     * ```
     */
    class touch_container {
      public:
        touch_container(base& p) : p_(p) {}
        inline Iterator begin() { return p_.begin(); }
        inline Iterator end() { return p_.end(); }

      private:
        base& p_;
    };

    /** Returns an iterable container for line `line` with this DIM. */
    touch_container operator()(line<D, T> line) {
        this->reset(line);
        return touch_container(*this);
    }

    /** Reset the DIM to a given line. */
    void reset(line<D, T> line) { reset_(line); }

    /** Obtain the current line of the DIM. */
    line<D, T> get_line() const { return line_; }

    /** Obtain a reference to the scanned volume. */
    volume<D> get_volume() const { return volume_; }

  protected:
    virtual ~base() = default;

    volume<D> volume_;
    line<D, T> line_;

  private:
    /* This idiom (keeping virtual functions private) is called NVI */
    virtual Iterator begin_() = 0;
    virtual Iterator end_() = 0;
    virtual void reset_(line<D, T> line) = 0;
};

} // namespace dim
} // namespace tomo
