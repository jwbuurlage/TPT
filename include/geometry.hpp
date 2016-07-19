#pragma once

#include <array>
#include <cassert>
#include <vector>

#include "math.hpp"

namespace tomo {

template <dimension Dimension, typename T = default_scalar_type>
struct line {
    using vecD = typename math::vec<Dimension, T>::type;

    line(vecD origin_in, vecD delta_in)
        : origin(origin_in), delta(delta_in) {}

    line(std::initializer_list<vecD> vectors)
        : line(*vectors.begin(), *(vectors.begin() + 1)) {
        // "trying to construct line with wrong number of vectors";
        assert(vectors.size() == 2);
    }

    // the direction vector
    vecD origin;

    // origin vector, first intersection of line with volume
    vecD delta;
};

// => Q: how do we deal with accelerators / distributed in this context
// Possible A: we still generate the begin/end points of the line on the
// CPU,
// but the 'kernel iterators' are implemented on the accelerator
template <dimension Dimension, typename T = default_scalar_type,
          class Iterator = void>
class geometry {
  public:
    geometry(int line_count) : line_count_(line_count) {}

    virtual Iterator begin() const;
    virtual Iterator end() const;

    int lines() const { return line_count_; }

    math::vec2<int> groups() const { return dimensions_; }

  protected:
    int line_count_;
    math::vec2<int> dimensions_;
};

} // namespace tomo
