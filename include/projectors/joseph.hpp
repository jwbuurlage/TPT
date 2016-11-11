#pragma once

namespace tomo {
namespace dim {

/** The iterator used for the Joseph DIM. */
template <typename T>
using joseph_iterator = typename std::vector<math::matrix_element<T>>::iterator;

/**
 * The Joseph DIM performs a single step along an axis, and the interpolates
 * between the other axes. A benefit of this technique is that there are no
 * 'shadowing non-zeros', i.e. non-zeros with the same indices.
 */
template <typename T>
class joseph: public base<2_D, T, joseph_iterator<T>> {
  public:
    /** Construct the DIM for a given volume. */
    joseph(volume<2_D> vol)
        : base<2_D, T, joseph_iterator<T>>(vol) {
        auto dims = this->volume_.dimensions();
        auto max_width = tomo::math::max_element<2_D, int>(dims);
        queue_.reserve((int)(2 * max_width));
    }

  private:
    void reset_(math::line<2_D, T> line) override {
        auto interpolate = [&](math::vec2<T> position, int axis) {
            math::vec2<int> closest_cell;
            closest_cell[1 - axis] = math::floor(position[1 - axis]);
            closest_cell[axis] =
                math::floor(position[axis] + 0.5 - math::epsilon<T>);

            math::vec2<int> cells[2] = {closest_cell, closest_cell};
            cells[1][axis] += 1;

            for (auto cell : cells) {
                auto distance = (T)1.0 - math::distance(cell[axis] + (T)0.5,
                                                        position[axis]);
                if (cell[axis] < this->volume_[axis] &&
                    distance > math::epsilon<T>) {
                    this->queue_.push_back(
                        {this->volume_.index({cell[0], cell[1]}), distance});
                }
            }
        };

        auto current_point = line.origin;
        int axis = (math::abs(line.delta.x) > math::abs(line.delta.y)) ? 1 : 0;

        auto step = line.delta / math::abs(line.delta[1 - axis]);
        current_point += (T)0.5 * step;

        while (math::inside<2_D, T>(current_point, this->volume_)) {
            interpolate(current_point, axis);
            current_point += step;
        }

        this->line_ = line;
    }

    void clear_() override {
        queue_.clear();
    }

    joseph_iterator<T> begin_() override { return queue_.begin(); }
    joseph_iterator<T> end_() override { return queue_.end(); }

    std::vector<math::matrix_element<T>> queue_;
};

} // namespace dim
} // namespace tomo
