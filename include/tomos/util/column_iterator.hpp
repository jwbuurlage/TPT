#pragma once

#include <array>
#include <limits>

#include "../geometry.hpp"
#include "../math/geometric.hpp"
#include "../projector.hpp"

namespace tomo {

template <dimension D, typename T>
class column {
  public:
    column(geometry::base<D, T>& geometry, dim::base<D, T>& kernel)
        : geometry_(geometry), kernel_(kernel), volume_(kernel_.get_volume()) {}

    const auto& operator()(math::vec<D, int> voxel) {
        using pix_iter = typename geometry::base<D, T>::pixel_iterator;

        values_.clear();

        // get corners of detector pixels
        auto corners = corners_(voxel);

        // project voxel to detector for each projection
        for (int p = 0; p < geometry_.projection_count(); ++p) {
            // if parallel/cone...
            auto shade = shadow_(corners, p);

            auto shape = math::vec<D - 1, int>(0);
            auto corner = geometry_.detector_corner(p);
            for (int d = 0; d < D - 1; ++d) {
                shape[d] = shade[d][1] - shade[d][0] + 1;
                corner += (T)shade[d][0] * geometry_.projection_delta(p)[d];
            }
            auto location = geometry_.source_location(p);
            auto delta = geometry_.projection_delta(p);
            auto parallel = geometry_.parallel();

            //            std::cout << "Projection: \n"
            //                      << "  s: " << math::to_string<D,
            //                      T>(location) << "\n"
            //                      << "  d: "
            //                      << math::to_string<D,
            //                      T>(geometry_.detector_corner(p))
            //                      << " ~ "
            //                      << math::to_string<D,
            //                      T>(geometry_.projection_delta(p)[0])
            //                      << "\n"
            //                      << "  v: " << math::to_string<D, int>(voxel)
            //                      << "\n"
            //                      << "shadow: [" << shade[0][0] << ", " <<
            //                      shade[0][1] << "]\n";

            auto detel_iterator =
                pix_iter(shape, corner, location, delta, parallel);

            // iterate over all these lines
            auto detels = math::product<D - 1, int>(shape);
            for (int i = 0; i < detels; ++i) {
                auto x = i;
                math::vec<D - 1, int> idx;
                for (int d = 0; d < D - 1; ++d) {
                    idx[d] = x % shape[d];
                    idx[d] += shade[d][0];
                    x /= shape[d];
                }
                int line_idx = idx[0];
                int offset = 1;
                for (int d = 1; d < D - 1; ++d) {
                    offset *= geometry_.projection_shape(p)[d];
                    line_idx += offset * idx[d];
                }
                line_idx += geometry_.offset(p);

                auto ray = *detel_iterator;
                // for each line, compute the matrix element
                // .. and add it to values_
                values_.push_back({line_idx, kernel_.matrix_value(ray, voxel)});
                ++detel_iterator;
            }
        }
        return *this;
    }

    auto begin() const { return values_.begin(); }
    auto end() const { return values_.end(); }

  private:
    auto corners_(math::vec<D, int> voxel) {
        const auto O = volume_.origin();
        const auto L = volume_.physical_lengths();
        const auto X = volume_.voxels();
        std::array<std::array<T, 2>, D> voxel_intervals;
        for (int d = 0; d < D; ++d) {
            voxel_intervals[d] = {O[d] + ((T)voxel[d] / X[d]) * L[d],
                                  O[d] + ((T)(voxel[d] + 1) / X[d]) * L[d]};
        }
        return math::cartesian_product<D, T>(voxel_intervals);
    }

    std::array<math::vec2<int>, D - 1>
    shadow_(std::array<math::vec<D, T>, math::pow(2, D)> corners, int p) {
        std::array<math::vec2<int>, D - 1> result;
        std::fill(result.begin(), result.end(),
                  math::vec2<int>{std::numeric_limits<int>::max(), -1});

        for (auto x : corners) {
            math::vec<D - 1, int> hyxel =
                cast_shadow_(x, geometry_.get_projection(p));
            for (int d = 0; d < D - 1; ++d) {
                result[d][0] = math::min(result[d][0], hyxel[d]);
                result[d][1] = math::max(result[d][1], hyxel[d]);
            }
        }
        for (int d = 0; d < D - 1; ++d) {
            result[d][0] = math::max(0, result[d][0] - 1);
            result[d][1] = math::min(geometry_.projection_shape(p)[d] - 1,
                                     result[d][1] + 1);
        }

        return result;
    }

    math::vec<D - 1, int> cast_shadow_(math::vec<D, T> x,
                                       geometry::projection<D, T> proj) {
        if (proj.parallel) {
            return math::parallel_project_hyxel(x, proj);
        } else {
            return math::project_hyxel(x, proj);
        }
    }

    geometry::base<D, T>& geometry_;
    dim::base<D, T>& kernel_;
    tomo::volume<D, T> volume_;

    std::vector<std::tuple<uint64_t, T>> values_;
};

} // namespace tomo
