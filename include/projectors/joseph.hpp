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
template <tomo::dimension D, typename T>
class joseph: public base<D, T, joseph_iterator<T>> {
  public:
    /** Construct the DIM for a given volume. */
    joseph(volume<D> vol)
        : base<D, T, joseph_iterator<T>>(vol) {
        auto dims = this->volume_.dimensions();
        auto max_width = tomo::math::max_element<D, int>(dims);
        queue_.reserve((int)(2 * max_width));
    }

  private:
    void reset_(math::line<D, T> line) override {
        auto interpolate = [&](math::vec<D, T> position, int axis) {
            math::vec<D, int> closest_cell;
            for (int d = 0; d < D; ++d) {
                if (d == axis)
                    continue;
                closest_cell[d] = math::floor(position[d]);
            }
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

        // the axis should correspond to the lowest delta component, for highest precision
        int axis = (math::abs(line.delta.x) > math::abs(line.delta.y)) ? 1 : 0;

        // FIXME actually implement this
        auto step = line.delta / math::abs(line.delta[1 - axis]);
        current_point += (T)0.5 * step;

        while (math::inside<D, T>(current_point, this->volume_)) {
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
