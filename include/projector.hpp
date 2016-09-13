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
    projector(volume<D> vol) : volume_(vol) {}

    Iterator begin() { return begin_(); }
    Iterator end() { return end_(); }

    /** This allows for code like:
      * for (auto l : g)
      *      for (auto touch : proj(l)) {
      *        // do something with touch things
      *      }
      */
    class touch_container {
      public:
        touch_container(projector& p) : p_(p) {}
        inline Iterator begin() { return p_.begin(); }
        inline Iterator end() { return p_.end(); }

      private:
        projector& p_;
    };

    touch_container operator()(line<D, T> line) {
        this->reset(line);
        return touch_container(*this);
    }

    void reset(line<D, T> line) { reset_(line); }

    line<D, T> get_line() const { return line_; }
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
