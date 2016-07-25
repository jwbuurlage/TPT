#pragma once

#include <cassert>

namespace tomo {

template <typename T>
using joseph_iterator = typename std::vector<math::matrix_element<T>>::iterator;

template <typename T>
class joseph_projector : public projector<2_D, T, joseph_iterator<T>> {
  public:
    joseph_projector(volume<2_D> vol)
        : projector<2_D, T, joseph_iterator<T>>(vol) {
        auto dims = this->volume_.dimensions();
        auto max_width = *std::max_element(dims.begin(), dims.end());
        queue_.reserve((int)(2 * max_width));
    }

    void reset(line<2_D, T> line) override {
        queue_.clear();

        auto interpolate = [&](math::vec2<T> position, int dimension) {
            math::vec2<int> closest_cell;
            closest_cell[1 - dimension] = math::floor(position[1 - dimension]);
            closest_cell[dimension] =
                math::floor(position[dimension] + 0.5 - math::epsilon);

            math::vec2<int> cells[2] = {closest_cell, closest_cell};
            cells[1][dimension] += 1;

            for (auto cell : cells) {
                auto distance = 1.0 - math::distance(cell[dimension] + 0.5,
                                                     position[dimension]);
                if (cell[dimension] < this->volume_[dimension] &&
                    distance > math::epsilon) {
                    this->queue_.push_back(
                        {this->volume_.index({cell[0], cell[1]}), distance});
                }
            }
        };

        auto current_point = line.origin;
        int dimension =
            (math::abs(line.delta.x) > math::abs(line.delta.y)) ? 1 : 0;

        auto step = line.delta / math::abs(line.delta[1 - dimension]);
        current_point += 0.5 * step;

        while (math::inside<2_D, T>(current_point, this->volume_)) {
            interpolate(current_point, dimension);
            current_point += step;
        }

        this->line_ = line;
    }

    joseph_iterator<T> begin() override { return queue_.begin(); }
    joseph_iterator<T> end() override { return queue_.end(); }

  private:
    std::vector<math::matrix_element<T>> queue_;
};

} // namespace tomo
