#pragma once

#include "tomo.hpp"

namespace tomo {

template <typename T>
using closest_iterator = typename std::vector<math::matrix_element<T>>::iterator;

template <typename T>
class closest_projector : public projector<2_D, T, closest_iterator<T>> {
  public:
    closest_projector(volume<2_D> vol)
        : projector<2_D, T, closest_iterator<T>>(vol) {
        auto dims = vol.dimensions();
        auto max_width = *std::max_element(dims.begin(), dims.end());
        queue_.reserve((int)(math::sqrt2 * max_width));
    }

    void reset(line<2_D, T> line) override {
        queue_.clear();

       auto current_point = line.origin;

        while (math::inside<2_D, T>(current_point, this->volume_)) {
            auto index = this->volume_.index(
                {(int)(current_point.x), (int)(current_point.y)});
            if (index >= 0 && index <= this->volume_.cells())
                queue_.push_back({index, (T)1.0});
            current_point += line.delta;
        }

        this->line_ = line;
    }

    closest_iterator<T> begin() override { return queue_.begin(); }
    closest_iterator<T> end() override { return queue_.end(); }

  private:
    std::vector<math::matrix_element<T>> queue_;
};

} // namespace tomo
