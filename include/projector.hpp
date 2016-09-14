#pragma once

#include <vector>

#include "common.hpp"
#include "geometry.hpp"
#include "logging.hpp"
#include "math.hpp"

namespace tomo {

/** A discrete integration method (DIM) for the ray inside a volume. */
template <dimension D, typename T, typename Iterator>
class projector {
  public:
    /** Construct a projector for a scanning volume. */
    projector(volume<D> vol) : volume_(vol) {}

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
        touch_container(projector& p) : p_(p) {}
        inline Iterator begin() { return p_.begin(); }
        inline Iterator end() { return p_.end(); }

      private:
        projector& p_;
    };

    /** Returns an iterable container for line `line` with this projector. */
    touch_container operator()(line<D, T> line) {
        this->reset(line);
        return touch_container(*this);
    }

    /** Reset the projector to a given line. */
    void reset(line<D, T> line) { reset_(line); }

    /** Obtain the current line of the projector. */
    line<D, T> get_line() const { return line_; }

    /** Obtain a reference to the scanned volume. */
    volume<D> get_volume() const { return volume_; }

  protected:
    virtual ~projector() = default;

    volume<D> volume_;
    line<D, T> line_;

  private:
    /* This idiom (keeping virtual functions private) is called NVI */
    virtual Iterator begin_() = 0;
    virtual Iterator end_() = 0;
    virtual void reset_(line<D, T> line) = 0;
};

} // namespace tomo
