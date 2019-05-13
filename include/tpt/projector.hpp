#pragma once

#include <type_traits>
#include <vector>

#include "common.hpp"
#include "geometry.hpp"
#include "logging.hpp"
#include "math.hpp"

namespace tpt {
namespace dim {

/**
 * A discrete integration method (DIM) for the ray inside a volume. Also called
 * a 'projector'
 */
template <dimension D, typename T>
class base {
  public:
    using problem_dimension = std::integral_constant<dimension, D>;
    using matrix_iterator =
        typename std::vector<math::matrix_element<T>>::iterator;

    /** Construct a projector for a scanning volume. */
    base(volume<D, T> vol) : volume_(vol) {}

    /** Obtain an iterator to the first voxel that is 'hit'. */
    matrix_iterator begin() { return begin_(); }

    /** Obtain an iterator beyond the final voxel that is 'hit'. */
    matrix_iterator end() { return end_(); }

    auto empty() { return queue_.empty(); }

    /**
     * This allows for code such as:
     * ```{.cpp}
     * for (auto l : g)
     *      for (auto element : proj(l)) {
     *        // do something with the element
     *      }
     * ```
     */
    auto& operator()(math::ray<D, T> incoming_ray) {
        // truncate to volume here first, then reset
        auto truncated_line = truncate_to_volume(incoming_ray, volume_);
        this->clear_();
        if (truncated_line) {
            auto line = truncated_line.value();

            this->reset_(line);
            this->line_ = line;
        }
        return *this;
    }

    /** Obtain the current line of the DIM. */
    math::line<D, T> get_line() const { return line_; }

    /** Obtain the scanned volume. */
    volume<D, T> get_volume() const { return volume_; }

    virtual T matrix_value(math::ray<D, T> ray, math::vec<D, int> voxel) = 0;

  protected:
    volume<D, T> volume_;
    math::line<D, T> line_;

    std::vector<math::matrix_element<T>> queue_;

  private:
    void clear_() { this->queue_.clear(); }
    matrix_iterator begin_() { return this->queue_.begin(); }
    matrix_iterator end_() { return this->queue_.end(); }

    /* This idiom (keeping virtual functions private) is called NVI */
    virtual void reset_(math::line<D, T> line) = 0;
};

} // namespace dim
} // namespace tpt
