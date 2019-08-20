#pragma once

#include <queue>
#include <random>
#include <vector>

#include "../../math/stringify.hpp"
#include "geometry.hpp"
#include "integration.hpp"

namespace tpt::grcb {

template <typename T>
std::tuple<cube<T>, cube<T>> split_at(cube<T> corners, T c, size_t d) {
    std::sort(corners.begin(), corners.end(),
              [&](auto& lhs, auto& rhs) { return lhs[d] < rhs[d]; });

    auto xs = corners;
    auto ys = corners;
    for (int i = 0; i < 4; ++i) {
        xs[4 + i][d] = c;
        ys[i][d] = c;
    }

    return {xs, ys};
}

/**
   input: arbitrarily ordered corners of a cube that is to be split, and the
   axis along which to split it midway. output: a pair containing the corners of
   the two parts, after splitting midway.

   Example in 2D:

   Before splitting:

   a - - - b
   |       |
   |       |
   |       |
   c - - - d

   After splitting (along x)

   a - e - b
   |   |   |
   |   |   |
   |   |   |
   c - f - d

   Function would return:
   (a, b, c, d) -> ((a, e, c, f), (e, b, f, d))
 */
template <typename T>
T split_midway(const geometry::base<3_D, T>& g [[maybe_unused]],
               cube<T> corners, size_t d) {
    std::sort(corners.begin(), corners.end(),
              [&](auto& lhs, auto& rhs) { return lhs[d] < rhs[d]; });

    auto midway =
        corners[0][d] + 0.5 * (corners[corners.size() - 1][d] - corners[0][d]);

    return midway;
}

template <typename T>
T split_equal_shadow(const geometry::base<3_D, T>& g, cube<T> corners,
                     size_t d) {
    std::sort(corners.begin(), corners.end(),
              [&](auto& lhs, auto& rhs) { return lhs[d] < rhs[d]; });
    auto xs = corners;
    auto ys = corners;

    auto min = 0.0;
    auto max = 1.0;
    while (max - min > 0.01f) {
        auto c = (max + min) * 0.5;
        auto midway = corners[0][d] +
                      c * (corners[corners.size() - 1][d] - corners[0][d]);
        xs = corners;
        ys = corners;
        for (int i = 0; i < 4; ++i) {
            xs[4 + i][d] = midway;
            ys[i][d] = midway;
        }

        // compute area of both shadows
        auto total_left = 0.0f;
        auto total_right = 0.0f;

        for (auto i = 0; i < g.projection_count(); ++i) {
            auto pi = g.get_projection(i);
            auto det = detector(pi);

            auto [pL, shadowL] = shadow(pi, xs);
            auto overlapL = std::vector<bg::model::polygon<vec2<T>>>{};
            bg::intersection(shadowL, det, overlapL);
            for (auto& sonic : overlapL) {
                total_left += bg::area(sonic);
            }

            auto [pR, shadowR] = shadow(pi, ys);
            auto overlapR = std::vector<bg::model::polygon<vec2<T>>>{};
            bg::intersection(shadowR, det, overlapR);
            for (auto& sonic : overlapR) {
                total_right += bg::area(sonic);
            }

            (void)pL;
            (void)pR;
        }

        if (math::abs(total_left - total_right) < math::epsilon<T>) {
            break;
        }
        if (total_left > total_right) {
            max = c;
        } else {
            min = c;
        }
    }

    auto c = (max + min) * 0.5;
    auto midway =
        corners[0][d] + c * (corners[corners.size() - 1][d] - corners[0][d]);

    return midway;
}

// template <typename T>
// T gradient(vec2<int> sample_axes,
//           std::vector<T> grid,
//           cube<T> v,
//           const geometry::base<3_D, T>& g,
//           const std::vector<math::vec3<T>>& sources,
//           T c,
//           int d,
//           T l,
//           T low,
//           T high)
//{
//    auto [vl, vr] = split_at(v, c, d);
//    auto left = load(vl, g);
//    auto right = load(vr, g);
//
//    auto sum = 0.0;
//    for (auto pnt : grid) {
//        sum += integrand(sources, {c, pnt}, {d, sample_axes[0]},
//        sample_axes[1], low, high);
//    }
//    auto avg = sum / grid.size();
//    auto fc = l * avg;
//
//    return 4.0 * fc * (left - right);
//}

// template <typename T>
// T split_point_gradient_descent(const cube<T> v, const geometry::base<3_D, T>&
// g, int d)
//{
//    auto sources = std::vector<math::vec3<T>>();
//    for (auto i = 0; i < g.projection_count(); ++i) {
//        sources.push_back(g.source_location(i));
//    }
//
//    auto [a, b] = min_max_cube(v);
//    auto c = (T)(a[d] + 0.5 * (b[d] - a[d]));
//    constexpr auto stepsize = (T)1e-4;
//
//    // we split at d, which leaves idx[0]
//    // for the numerical integration, and idx[1] for the analytiacl integral
//    auto idx = math::restrict<3_D, T>({0, 1, 2}, d);
//    auto grid = std::vector<T>();
//    auto k = 50;
//    auto l = (b[idx[0]] - a[idx[0]]);
//    auto h = l / (T)(k - 1);
//    for (auto i = 0; i < k; ++i) {
//        grid.push_back(a[idx[0]] + i * h);
//    }
//
//    for (int step = 0; step < 10; ++step) {
//        auto delta = gradient(idx, grid, v, g, sources, c, d, l, a[idx[1]],
//        b[idx[1]]);
//
//        c -= stepsize * delta;
//        if (std::fabs(delta) < 1e-3) {
//            break;
//        }
//
//        // box constraint on c
//        c = std::max(c, a[d]);
//        c = std::min(c, b[d]);
//    }
//
//    return c;
//}

template <typename T>
T split_point_binary_search(precision level, method kind,
                            const geometry::base<3_D, T>& g, cube<T> v,
                            size_t d) {
    auto sources = std::vector<math::vec3<T>>();
    for (auto i = 0; i < g.projection_count(); ++i) {
        sources.push_back(g.source_location(i));
    }

    auto [a, b] = min_max_cube<T>(v);
    auto c = (T)(a[d] + 0.5 * (b[d] - a[d]));

    auto min = a[d];
    auto max = b[d];
    const auto tol = 1e-8;

    // max 100 steps
    for (int step = 0; step < 100; ++step) {
        auto [vl, vr] = split_at(v, c, d);
        auto left = load(vl, g, level, kind);
        auto right = load(vr, g, level, kind);

        auto delta = right - left;

        if ((std::fabs((2 * delta) / (left + right)) < tol) ||
            (((max - min) / (b[d] - a[d])) < tol)) {
            break;
        }

        if (delta > 0) {
            min = c;
        } else {
            max = c;
        }
        c = 0.5f * (max + min);
    }

    return c;
}

// SMART: Split at Median Aggregate for Ray Tracing
template <typename T>
T split_smart(precision level, const geometry::base<3_D, T>& g, cube<T> v,
              size_t d) {
    auto sources = std::vector<math::vec3<T>>();
    for (auto i = 0; i < g.projection_count(); ++i) {
        sources.push_back(g.source_location(i));
    }

    auto [a, b] = min_max_cube<T>(v);

    auto rho = 1e-8;
    auto n_rho = 10;

    switch (level) {
    case precision::low: {
        rho = 1e-3;
        n_rho = 10;
        break;
    }
    case precision::medium: {
        rho = 1e-5;
        n_rho = 50;
        break;
    }
    case precision::high: {
        rho = 1e-8;
        n_rho = 100;
        break;
    }
    default:
        break;
    }

    auto rd = std::random_device();
    auto engine = std::mt19937(rd());
    auto distribution = std::uniform_real_distribution<T>((T)0, (T)1);

      auto f = [&](auto x) {
        auto sum = (T)0;
        for (auto i = 0; i < g.projection_count(); ++i) {
            auto s = g.source_location(i);
            auto pi = g.get_projection(i);
            auto proj_point = world_to_detector(pi, project(pi, x).value());
            auto hds = (T)0.5 * pi.detector_size;
            bool outside = false;
            for (int d = 0; d < 2; ++d) {
                if (proj_point[d] < (T)-hds[d] || proj_point[d] > (T)hds[d]) {
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

    using point = std::pair<T, T>;

    auto sample = [&, a, b, d]() {
        auto x = vec3<T>();
        for (auto k = 0u; k < 3; ++k) {
            x[k] = a[k] + (b[k] - a[k]) * distribution(engine);
        }
        return point{x[d], f(x)};
    };

    // 1. bootstrap the process, take two random samples, and find the C
    auto l = sample();
    auto r = sample();
    if (std::get<0>(l) > std::get<0>(r)) {
        std::swap(l, r);
    }

    auto left =
        std::priority_queue<point, std::vector<point>, std::less<point>>();
    left.push(l);
    auto right =
        std::priority_queue<point, std::vector<point>, std::greater<point>>();
    right.push(r);

    // left is:
    // l * s * (c - a_d)
    // right is
    // r * s * (b_d - c)
    // so to find equal c:
    // l * (c - ad) = r * (b_d - c)
    // l * c - l * a_d = r * b_d - r * c
    // (l + r) * c = r * b_d + l * a_ad
    // c = ( r * b_d + l * a_d ) / (l + r)

    // iteratively update c
    auto sum_l = std::get<1>(l);
    auto sum_r = std::get<1>(r);

    auto compute_c = [a, b, d](auto max_l, auto min_r, auto avg_l, auto avg_r) {
        auto c = (avg_r * b[d] + avg_l * a[d]) / (avg_l + avg_r);

        // should be boxed between l and r
        if (c > min_r) {
            c = min_r;
        }
        if (c < max_l) {
            c = max_l;
        }
        return c;
    };

    auto integrals = [a, b, d](auto c, auto avg_l, auto avg_r) {
        auto L = avg_l * (c - a[d]);
        auto R = avg_r * (b[d] - c);

        return std::pair<T, T>{L, R};
    };

    auto c = compute_c(std::get<0>(l), std::get<0>(r), sum_l, sum_r);

    auto k_rho = 0;

    auto cs = std::queue<T>();

    auto max_iters = 10000;
    auto iter = 0;
    while (true) {
        // 1) generate a new sample
        auto [x, fx] = sample();
        auto max_l = std::get<0>(left.top());
        auto min_r = std::get<0>(right.top());
        // 2) see if it gets inserted left or right, update sum of that part
        if (x < max_l) {
            left.push(point{x, fx});
            sum_l += fx;
            max_l = std::get<0>(left.top());
        } else {
            right.push(point{x, fx});
            sum_r += fx;
            min_r = std::get<0>(right.top());
        }

        // swapping phase
        // 3) compare left and right parts, (with computed c) two cases e.g.
        // left > right:
        auto cc =
            compute_c(max_l, min_r, sum_l / left.size(), sum_r / right.size());
        auto [L, R] = integrals(cc, sum_l / left.size(), sum_r / right.size());
        //        std::cout << "Iteration: c = " << cc << ", L = " << L << ", R
        //        = " << R
        //                  << ". left.size() = " << left.size()
        //                  << ", right.size() = " << right.size() << "\n";

        if (L > R) {
            while (true) {
                // always let one remain in L
                if (left.size() == 1) {
                    break;
                }
                // virtually move top point from left to right
                auto [x, fx] = left.top();

                // peek underneath current top
                left.pop();
                auto max_l_new = std::get<0>(left.top());
                left.push({x, fx});

                auto min_r_new = x;
                auto avg_l_new = (sum_l - fx) / (left.size() - 1);
                auto avg_r_new = (sum_r + fx) / (right.size() + 1);
                auto c_new =
                    compute_c(max_l_new, min_r_new, avg_l_new, avg_r_new);
                auto [L_new, R_new] = integrals(c_new, avg_l_new, avg_r_new);
                if (std::fabs(R_new - L_new) < std::fabs(R - L)) {
                    // accept and try again
                    left.pop();
                    right.push({x, fx});
                    sum_l -= fx;
                    sum_r += fx;
                    L = L_new;
                    R = R_new;
                } else {
                    break;
                }
            }
        } else {
            while (true) {
                // always let one remain in R
                if (right.size() == 1) {
                    break;
                }
                // virtually move top point from right to left
                auto [x, fx] = right.top();
                auto max_l_new = x;

                // peek underneath current top
                right.pop();
                auto min_r_new = std::get<0>(right.top());
                right.push({x, fx});

                auto avg_l_new = (sum_l + fx) / (left.size() + 1);
                auto avg_r_new = (sum_r - fx) / (right.size() - 1);

                auto c_new =
                    compute_c(max_l_new, min_r_new, avg_l_new, avg_r_new);
                auto [L_new, R_new] = integrals(c_new, avg_l_new, avg_r_new);

                if (std::fabs(R_new - L_new) < std::fabs(R - L)) {
                    // accept and try again
                    right.pop();
                    left.push({x, fx});
                    sum_l += fx;
                    sum_r -= fx;
                    L = L_new;
                    R = R_new;
                } else {
                    break;
                }
            }
        }

        // find new C
        max_l = std::get<0>(left.top());
        min_r = std::get<0>(right.top());
        auto c_new =
            compute_c(max_l, min_r, sum_l / left.size(), sum_r / right.size());

        cs.push(c_new);

        // did our c change compared to n_rho ago?
        if (cs.size() > (size_t)n_rho) {
            if (std::fabs(c_new - cs.front()) < rho) {
                k_rho++;
            } else {
                k_rho = 0;
            }
            cs.pop();
        }

        c = c_new;

        // if not for N_rho times, break here
        if (k_rho >= n_rho) {
            break;
        }

        ++iter;
        if (iter > max_iters) {
            break;
        }
    }

    return c;
}

template <typename T>
T split_smart_simple(precision level, const geometry::base<3_D, T>& g,
                     cube<T> v, size_t d) {

    auto n = std::map<precision, int>{{precision::low, 2500},
                                      {precision::medium, 5000},
                                      {precision::high, 20000}}[level];

    auto sources = std::vector<math::vec3<T>>();
    for (auto i = 0; i < g.projection_count(); ++i) {
        sources.push_back(g.source_location(i));
    }

    auto [a, b] = min_max_cube<T>(v);

    auto rd = std::random_device();
    auto engine = std::mt19937(rd());
    auto distribution = std::uniform_real_distribution<T>((T)0, (T)1);

    auto f = [&](auto x) {
        auto sum = (T)0;
        for (auto i = 0; i < g.projection_count(); ++i) {
            auto s = g.source_location(i);
            auto pi = g.get_projection(i);
            auto proj_point = world_to_detector(pi, project(pi, x).value());
            auto hds = (T)0.5 * pi.detector_size;
            bool outside = false;
            for (int d = 0; d < 2; ++d) {
                if (proj_point[d] < (T)-hds[d] || proj_point[d] > (T)hds[d]) {
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

    using point = std::pair<T, T>;

    auto sample = [&, a, b, d]() {
        auto x = vec3<T>();
        for (auto k = 0u; k < 3; ++k) {
            x[k] = a[k] + (b[k] - a[k]) * distribution(engine);
        }
        return point{x[d], f(x)};
    };

    auto compute_c = [a, b, d](auto max_l, auto min_r, auto avg_l, auto avg_r) {
        auto c = (avg_r * b[d] + avg_l * a[d]) / (avg_l + avg_r);

        // should be boxed between l and r
        if (c > min_r) {
            c = min_r;
        }
        if (c < max_l) {
            c = max_l;
        }
        return c;
    };

    auto xs = std::vector<point>();
    for (int i = 0; i < n; ++i) {
        xs.push_back(sample());
    }

    std::sort(xs.begin(), xs.end());

    auto sum = (T)0;
    for (auto [x, fx] : xs) {
        sum += fx;
    }

    auto sum_l = (T)0;
    auto r = (T)0;
    int k = 0;
    for (; k < n - 1; ++k) {
        auto [x, fx] = xs[k];
        auto [y, fy] = xs[k + 1];
        // unused
        (void)x, (void)fy;

        // check turning point where L > R
        sum_l += fx;
        r = y;

        auto L = (sum_l / (k + 1)) * (r - a[d]);
        auto R = ((sum - sum_l) / (n - (k + 1))) * (b[d] - r);
        if (L > R) {
            break;
        }
    }

    auto c = compute_c(std::get<0>(xs[k]), std::get<0>(xs[k + 1]),
                       sum_l / (k + 1), ((sum - sum_l) / (n - (k + 1))));

    return c;
}

} // namespace tpt::grcb
