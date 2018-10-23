#pragma once

#include "../math.hpp"
#include "../projector.hpp"

namespace tomo {
namespace dim {

/**
 * The Siddon DIM is an efficient implementation of a line kernel.
 */
template <typename T>
class siddon : public base<3_D, T> {
    static constexpr auto D = 3_D;

  public:
    /** Construct the DIM for a given volume. */
    siddon(volume<D, T> vol) : base<D, T>(vol) {}

    T matrix_value(math::ray<D, T> ray, math::vec<D, int> voxel) override {
        (void)ray;
        (void)voxel;
        // TODO implement
        return -1;
    }

  private:
    void reset_(math::line<D, T> line) override {
        auto N = math::vec<D, T>(this->volume_.voxels());

        // compute entry/exit points
        auto start_end = math::aabb_intersection<D, T>(
            line.origin - line.delta, line.origin + N[0] * 2.0f * line.delta,
            this->volume_.physical_lengths(), this->volume_.origin());
        if (!start_end) {
            return;
        }
        auto [p1, p2] = start_end.value();

        // convert to voxel space
        p1 = to_voxel_(p1);
        p2 = to_voxel_(p2);

        // FIXME always one of the dimension constant..
        // What is the proper way to fix this
        for (auto d = 0u; d < D; ++d) {
          if (math::approx_equal(p1[d], p2[d])) {
            p1[d] = 0;
            p2[d] = 1;
          }
        }

        auto alpha = [p1, p2](auto i) { return (i - p1) / (p2 - p1); };

        auto min = alpha(math::vec<D, T>{0, 0, 0});
        auto max = alpha(N - math::vec<D, T>{1, 1, 1});

        // swap alphas with min/max
        auto swapped_min = math::min(min, max);
        auto swapped_max = math::max(min, max);

        auto alpha_min = math::min_element<D, T>(swapped_min);
        auto alpha_max = math::max_element<D, T>(swapped_max);

        std::cout << "1. " << glm::to_string(p1) << " -> " << glm::to_string(p2)
                  << "\n";
        std::cout << "2. " << glm::to_string(min) << " -> "
                  << glm::to_string(max) << "\n";
        std::cout << "3. " << glm::to_string(swapped_min) << " -> "
                  << glm::to_string(swapped_max) << "\n";
        std::cout << "4. " << alpha_min << " -> " << alpha_max << "\n";
        std::cout << "---\n";

        // std::vector<math::matrix_element<T>> queue_; {index, value}
        // TODO implement
        this->queue_.push_back({0, (T)1});
    }

    auto to_voxel_(math::vec<D, T> x) {
        return (x - this->volume_.origin()) / this->volume_.physical_lengths() *
               math::vec<D, T>(this->volume_.voxels());
    }
};

} // namespace dim
} // namespace tomo
