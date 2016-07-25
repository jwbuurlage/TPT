#pragma once

#include "../common.hpp"
#include "../geometry.hpp"

namespace tomo {

template <dimension D, typename T = default_scalar_type>
using list_line_iterator = typename std::vector<line<D, T>>::const_iterator;

template <dimension D, typename T = default_scalar_type>
class list_geometry : public geometry<D, T, list_line_iterator<D, T>> {
  public:
    using iterator_type = list_line_iterator<D, T>;

    list_geometry(std::vector<line<D, T>>& lines)
        : geometry<D, T, iterator_type>(lines.size()), lines_(lines) {}

    list_geometry(std::vector<line<D, T>>&& lines)
        : geometry<D, T, iterator_type>(lines.size()), lines_(lines) {}

    iterator_type begin() const override { return lines_.begin(); }
    iterator_type end() const override { return lines_.end(); }

    void set_dimensions(math::vec2<T> dimensions) {
        this->dimensions_ = dimensions;
    }

  private:
    std::vector<line<D, T>> lines_;
};

} // namespace tomo
