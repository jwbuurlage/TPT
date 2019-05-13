#pragma once

#include <cmath>
#include <random>
#include <vector>

#include "datatypes.hpp"
#include "../../math.hpp"
#include "../../geometry.hpp"

namespace tpt::grcb {

enum class method { two_dim, three_dim, cone };
enum class precision { low, medium, high };

template <uint32_t D, typename T>
class integrator {
  public:
    integrator(precision level) : level_(level) {}

    virtual T integrate(box<D, T> domain, std::function<T(vec<D, T>)> f) = 0;

  protected:
    precision level_;
};

template <uint32_t D, typename T>
class simple_sampling : public integrator<D, T> {
  public:
    simple_sampling(precision level) : integrator<D, T>(level) {
        auto n = 0u;
        switch (this->level_) {
        case precision::low: {
            n = 100u;
            break;
        }
        case precision::medium: {
            n = 1000u;
            break;
        }
        case precision::high: {
            n = 10000u;
            break;
        }
        };

        auto rd = std::random_device();
        auto engine = std::mt19937(rd());
        auto distribution = std::uniform_real_distribution<T>((T)0, (T)1);

        grid_.resize(n);
        for (auto i = 0u; i < n; ++i) {
            for (auto d = 0u; d < D; ++d) {
                grid_[i][d] = distribution(engine);
            }
        }
    }

    virtual T integrate(box<D, T> domain,
                        std::function<T(vec<D, T>)> f) override {
        auto sum = (T)0;
        for (auto u : grid_) {
            sum += f(domain[0] + (domain[1] - domain[0]) * u);
        }

        auto area = math::product<D, T>(domain[1] - domain[0]);
        return area * (sum / grid_.size());
    }

  private:
    std::vector<vec<D, T>> grid_;
};

template <typename T>
T integrand(const std::vector<vec3<T>>& sources, vec2<T> sample_point,
            vec2<int> sample_axes, int d, T low, T high) {
    auto x = sample_point[0];
    auto y = sample_point[1];

    auto sum = 0.0;
    for (auto s : sources) {
        auto a = std::sqrt((x - s[sample_axes[0]]) * (x - s[sample_axes[0]]) +
                           (y - s[sample_axes[1]]) * (y - s[sample_axes[1]]));
        sum += (1.0 / a) *
               (std::atan2((high - s[d]), a) - std::atan2((low - s[d]), a));
    }

    return sum;
}

// TODO how to pass integration strategy, precsision
template <typename T>
T load(cube<T> v, const geometry::base<3_D, T>& g, precision level,
       method kind) {
    auto sources = std::vector<math::vec3<T>>();
    for (auto i = 0; i < g.projection_count(); ++i) {
        sources.push_back(g.source_location(i));
    }

    switch (kind) {
    case method::two_dim: {
        auto ninter = simple_sampling<2_D, T>(level);
        auto [a, b] = min_max_cube<T>(v);
        return ninter.integrate(
            box<2_D, T>{vec2<T>{a[0], a[1]}, vec2<T>{b[0], b[1]}}, [&](auto x) {
                return integrand(sources, x, {0, 1}, 2, a[2], b[2]);
            });
        break;
    }

    case method::three_dim: {
        auto ninter = simple_sampling<3_D, T>(level);
        auto [a, b] = min_max_cube<T>(v);
        return ninter.integrate({a, b}, [&](auto x) {
            auto sum = (T)0;
            for (auto s : sources) {
                sum += (T)1 / math::sum<3_D, T>((x - s) * (x - s));
            }
            return sum;
        });
        break;
    }

    case method::cone: {
        auto ninter = simple_sampling<3_D, T>(level);
        auto [a, b] = min_max_cube<T>(v);
        return ninter.integrate({a, b}, [&](auto x) {
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
        });
        break;
    }

    default: {
        std::cout << "Warning: Unknown integration method called.\n";
        break;
    }
    }

    return (T)0;
}

} // namespace tpt::grcb
