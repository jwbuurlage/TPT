#pragma once

#include "../common.hpp"
#include "../geometry.hpp"

namespace tomo {
namespace geometry {

/**
 * Geometry defined by a list of rays.
 *
 * \tparam D the dimension of the volume.
 * \tparam T the scalar type to use
 */
template <dimension D, typename T = default_scalar_type>
class list : public base<D, T, list<D, T>> {
  public:
    /** Construct a list geometry by copying a vector of lines. */
    list(const std::vector<math::ray<D, T>>& lines)
        : base<D, T, list<D, T>>(lines.size()), lines_(lines) {}

    /** Construct a list geometry by moving a vector of lines. */
    list(std::vector<math::ray<D, T>>&& lines)
        : base<D, T, list<D, T>>(lines.size()), lines_(lines) {}

    /** Obtain the i-th line of the geometry. */
    inline math::ray<D, T> get_line(int i) const { return lines_[i]; }

    /**
     * Set the dimensions of the list geometry.
     * \see `tomo::geometry::groups()`
     */
    // FIXME this should possibly be renamed to `set_groups`.
    void set_dimensions(math::vec2<T> dimensions) {
        this->dimensions_ = dimensions;
    }

  private:
    std::vector<math::ray<D, T>> lines_;
};

} // namespace geometry
} // namespace tomo
