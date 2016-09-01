#pragma once

#include "tomo.hpp"

namespace tomo {

template <typename T>
using closest_iterator =
    typename std::vector<math::matrix_element<T>>::iterator;

template <dimension D, typename T>
class closest_projector : public projector<D, T, closest_iterator<T>> {
  public:
    closest_projector(volume<D> vol)
        : projector<D, T, closest_iterator<T>>(vol) {
        auto dims = vol.dimensions();
        auto max_width = *std::max_element(dims.begin(), dims.end());
        queue_.reserve((int)(math::sqrt(D) * max_width));
    }

    void reset(line<D, T> line) override {
        queue_.clear();

        auto current_point = line.origin;

        while (math::inside<D, T>(current_point, this->volume_)) {
			std::array<int, D> xs;
			for (int i = 0; i < D; ++i) {
				xs[i] = (int)current_point[i];
			}
            auto index = this->volume_.index(xs);
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
