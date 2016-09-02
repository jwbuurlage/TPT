#pragma once

#include <algorithm>

#include "../projector.hpp"
#include "../volume.hpp"

namespace tomo {

// bi- or trilinear projector
template <dimension D, typename T, typename Iterator>
class dlinear_projector : public projector<D, T, Iterator> {
  public:
    dlinear_projector(volume<D> vol) : projector<D, T, Iterator>(vol) {
        auto dims = this->volume_.dimensions();
        auto max_width = *std::max_element(dims.begin(), dims.end());
        queue_.reserve((int)(math::sqrt2<T> * 4 * max_width));
    }

    void reset(line<D, T> line) override {
        queue_.clear();
        auto current_point = line.origin;

        while (math::inside<D, T>(current_point, this->volume_)) {
            math::interpolate<D, T>(current_point, this->volume_, queue_);
            current_point += line.delta;
        }

        this->line_ = line;
    }

    Iterator begin() override { return queue_.begin(); }
    Iterator end() override { return queue_.end(); }

  private:
    std::vector<math::matrix_element<T>> queue_;
};

template <dimension D, typename T = default_scalar_type>
using linear_projector =
    dlinear_projector<D, T,
                      typename std::vector<math::matrix_element<T>>::iterator>;

} // namespace tomo
