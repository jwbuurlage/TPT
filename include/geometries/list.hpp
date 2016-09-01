#pragma once

#include "../common.hpp"
#include "../geometry.hpp"

namespace tomo {

template <dimension D, typename T = default_scalar_type>
class list_geometry : public geometry<D, T, list_geometry<D, T>> {
  public:
    list_geometry(std::vector<line<D, T>>& lines)
        : geometry<D, T, list_geometry<D, T>>(lines.size()), lines_(lines) {}

    list_geometry(std::vector<line<D, T>>&& lines)
        : geometry<D, T, list_geometry<D, T>>(lines.size()), lines_(lines) {}

    inline line<D, T> get_line(int i) const { return lines_[i]; }

    void set_dimensions(math::vec2<T> dimensions) {
        this->dimensions_ = dimensions;
    }

  private:
    std::vector<line<D, T>> lines_;
};

} // namespace tomo
