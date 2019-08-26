#pragma once

#include <limits>
#include <memory>

#include "geometry.hpp"

namespace tpt::grcb {

/** A tree. */
template <typename T>
struct node {
    std::unique_ptr<node<T>> left;
    int a;
    T c;
    std::unique_ptr<node<T>> right;
};

/** Output a tree. */
template <typename T>
void print_tree(const std::unique_ptr<node<T>>& n, std::ostream& out) {
    if (!n) {
        out << "[]";
        return;
    }
    out << "[[" << n->a << ", " << n->c << "], ";
    print_tree(n->left, out);
    out << ", ";
    print_tree(n->right, out);
    out << "]";
}

/**
 * Note: the unit cube (origin, lengths) is for outputting the tree in neutral
 * coordinates, it is not the bounding box of the volume.
 */

template <typename T>
std::unique_ptr<node<T>> split(math::vec3<T> origin, math::vec3<T> lengths,
                               cube<T> v, const geometry::base<3_D, T>& g,
                               int q, bool use_gradient, T mag,
                               split_function<T> splitter, int depth = 0) {
    auto abar = -1;
    auto Abar = std::numeric_limits<T>::max();
    auto cbar = (T)0;

    for (auto a = 0u; a < 3; ++a) {
        auto c = splitter(g, v, a);

        auto [vL, vR] = split_at(v, c, a);

        auto A = 0.0;
        for (auto i = 0; i < g.projection_count(); ++i) {
            A += overlap(g.get_projection(i), vL, vR);
        }

        if (use_gradient) {
            A += mag * cross_section_area<T>(v, a);
        }

        if (A < Abar) {
            abar = a;
            Abar = A;
            cbar = c;
        }
    }

    // update the 'neutral cube' output
    // std::cout << "cbar: " << cbar << "\n";
    auto [a, b] = min_max_cube<T>(v);
    auto neutral_split = (cbar - a[abar]) / (b[abar] - a[abar]);

    auto lengthsL = lengths;
    auto lengthsR = lengths;
    lengthsL[abar] = neutral_split * lengths[abar];
    lengthsR[abar] = ((T)1 - neutral_split) * lengths[abar];
    auto originR = origin;
    originR[abar] += lengthsL[abar];

    if (depth + 1 == q) {
        return std::unique_ptr<node<T>>(
            new node<T>{nullptr, abar, originR[abar], nullptr});
    } else {
        auto [vL, vR] = split_at(v, cbar, abar);
        return std::unique_ptr<node<T>>(
            new node<T>{split(origin, lengthsL, vL, g, q, use_gradient, mag,
                              splitter, depth + 1),
                        abar, originR[abar],
                        split(originR, lengthsR, vR, g, q, use_gradient, mag,
                              splitter, depth + 1)});
    }
}

/** Perform GRCB with midway splits */
template <typename T>
std::unique_ptr<node<T>> partition(volume<3_D, T> v,
                                   const geometry::base<3_D, T>& g, int q,
                                   std::vector<flag> flags) {
    auto vcorners = corners(v);

    auto det_size = g.get_projection(0).detector_size;
    auto mag = cross_section_area<T>(vcorners, 0) / (det_size[0] * det_size[1]);

    auto passed = [&](auto the_flag) -> bool {
        auto pos = std::find(std::begin(flags), std::end(flags), the_flag);
        return pos != flags.end();
    };

    auto use_gradient = passed(flag::gradient_volume);
    auto split_midway_passed [[maybe_unused]] = passed(flag::split_midway);
    auto split_equal_shadow_passed = passed(flag::split_equal_shadow);
    auto split_equal_load_passed = passed(flag::split_equal_load);
    auto split_smart_passed = passed(flag::split_smart);
    auto split_smart_simple_passed = passed(flag::split_smart_simple);

    auto kind = method::two_dim;
    if (passed(flag::integration_three_dim)) {
        kind = method::three_dim;
    } else if (passed(flag::integration_cone)) {
        kind = method::cone;
    }
    auto level = precision::low;
    if (passed(flag::precision_medium)) {
        level = precision::medium;
    } else if (passed(flag::precision_high)) {
        level = precision::high;
    }

    auto splitter = split_function<T>();
    if (split_equal_load_passed) {
        splitter = [=](auto& g, auto v, auto d) {
            return split_point_binary_search<T>(level, kind, g, v, d);
        };
    } else if (split_equal_shadow_passed) {
        splitter = split_equal_shadow<T>;
    } else if (split_smart_passed) {
        splitter = [=](auto& g, auto v, auto d) {
            return split_smart<T>(level, g, v, d);
        };
    } else if (split_smart_simple_passed) {
        static auto ns =
            std::map<precision, int>{{precision::low, 100'000},
                                     {precision::medium, 1'000'000},
                                     {precision::high, 10'000'000}};
        auto n = ns[level];

        auto sources = std::vector<math::vec3<T>>(g.projection_count());
        for (auto i = 0; i < g.projection_count(); ++i) {
            sources[i] = g.source_location(i);
        }

        auto [a, b] = min_max_cube<T>(vcorners);

        auto f = [&](auto x) {
            auto sum = (T)0;
            for (auto i = 0; i < g.projection_count(); ++i) {
                auto s = g.source_location(i);
                auto pi = g.get_projection(i);
                auto proj_point = world_to_detector(pi, project(pi, x).value());
                auto hds = (T)0.5 * pi.detector_size;
                bool outside = false;
                for (int d = 0; d < 2; ++d) {
                    if (proj_point[d] < (T)-hds[d] ||
                        proj_point[d] > (T)hds[d]) {
                        outside = true;
                        break;
                    }
                }
                if (!outside) {
                    sum += (T)1 / math::sum<3_D, T>((x - s) * (x - s));
                }
            }
            return sum;
        };

        using sample = std::pair<math::vec3<T>, T>;

        auto rd = std::random_device();
        auto engine = std::mt19937(rd());
        auto distribution = std::uniform_real_distribution<T>((T)0, (T)1);

        auto samples = std::vector<sample>(n);
        auto idx = 0u;
        for (auto i = 0; i < n; ++i) {
            auto x = math::vec3<T>{};
            for (auto l = 0; l < 3; ++l) {
                x[l] = a[l] + (b[l] - a[l]) * distribution(engine);
            }
            samples[idx++] = sample{x, f(x)};
        }

        // preprocess grid points and pass to split smart simple
        splitter = [=](auto&, auto v, auto d) {
            return split_smart_simple<T>(v, d, samples);
        };
    } else {
        splitter = split_midway<T>;
    }

    // Note: the unit cube is for outputting the tree in neutral coordinates, it
    // is not the bounding box of the volume
    return split(math::vec3<T>{0, 0, 0}, math::vec3<T>{1, 1, 1}, vcorners, g, q,
                 use_gradient, mag, splitter);
}

} // namespace tpt::grcb
