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
    projector(volume<D> vol) : volume_(vol), line_({{0, 0}, {0, 0}}) {}

    virtual Iterator begin() = 0;
    virtual Iterator end() = 0;

    virtual void reset(line<D, T> line)  = 0;

    line<D, T> get_line() const { return line_; }
    volume<D> get_volume() const { return volume_; }

  protected:
    volume<D> volume_;
    line<D, T> line_;
};

} // namespace tomo
