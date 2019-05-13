#pragma once

#include "../geometry.hpp"

namespace tpt {
namespace geometry {

/**
 * A geometry
 *
 * \tparam D the dimension of the volume.
 * \tparam T the scalar type to use
 */
template <typename T>
class custom : public base<3_D, T> {
  public:
    custom(std::vector<projection<3_D, T>> projection_list)
        : base<3_D, T>(projection_list.size()),
          projection_list_(projection_list) {}

    math::vec<3_D - 1, int> projection_shape(int i) const override {
        return projection_list_[i].detector_shape;
    }

    math::vec<3_D, T> detector_corner(int i) const override {
        math::vec<3_D, T> result = projection_list_[i].detector_location;
        for (int d = 0; d < 3_D - 1; ++d) {
            result -= (T)0.5 * (projection_list_[i].detector_size[d] *
                                projection_list_[i].detector_tilt[d]);
        }
        return result;
    }

    math::vec<3_D, T> source_location(int i) const override {
        return projection_list_[i].source_location;
    }

    projection<3_D, T> get_projection(int idx) const override {
        return projection_list_[idx];
    }

    std::array<math::vec<3_D, T>, 3_D - 1>
    projection_delta(int i) const override {
        auto axes = projection_list_[i].detector_tilt;
        std::array<math::vec<3_D, T>, 3_D - 1> result;

        for (int d = 0; d < 3_D - 1; ++d) {
            result[d] = (projection_list_[i].detector_size[d] /
                         projection_list_[i].detector_shape[d]) *
                        axes[d];
        }

        return result;
    }

  private:
    std::vector<projection<3_D, T>> projection_list_;
};

} // namespace geometry
} // namespace tpt
