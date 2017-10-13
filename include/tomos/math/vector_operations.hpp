#pragma once

#include "constants.hpp"
#include "vector.hpp"

#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>

namespace tomo {

using namespace tomo::literals;

namespace math {

/** Compute the absolute value of an object of type T. */
template <typename T>
auto abs(T obj) {
    return glm::abs(obj);
}

/** Check if two real-valued numbers are approximately equal. */
template <typename T, typename V>
bool approx_equal(T lhs, T rhs, V max_rel_diff, V max_abs_diff) {
    return abs((lhs - rhs) / rhs) <= max_rel_diff ||
           (max_abs_diff != 0 && abs(lhs - rhs) <= max_abs_diff);
}

/** Check if two real-valued numbers are approximately equal. */
template <typename T>
bool approx_equal(T lhs, T rhs) {
    return approx_equal(lhs, rhs, (T)1e-2, epsilon<T>);
}

/** Compute the cosine for an object of type T. */
template <typename T>
auto cos(T obj) {
    return glm::cos(obj);
}

/** Compute the max for an object of type T. */
template <typename T>
auto max(T a, T b) {
    return glm::max(a, b);
}

/** Compute the max for an object of type T. */
template <typename T>
auto min(T a, T b) {
    return glm::min(a, b);
}

/** Compute the sine for an object of type T. */
template <typename T>
auto sin(T obj) {
    return glm::sin(obj);
}

/** Compute the arcsine for an object of type T. */
template <typename T>
auto asin(T obj) {
    return glm::asin(obj);
}

/** Compute the exponential of an object of type T. */
template <typename T>
auto exp(T obj) {
    return glm::exp(obj);
}

/** Floor an object of type T. */
template <typename T>
auto floor(T obj) {
    return glm::floor(obj);
}

/** Ceil an object of type T. */
template <typename T>
auto ceil(T obj) {
    return glm::ceil(obj);
}

/** Take the logarithm of an object of type T. */
template <typename T>
auto log(T obj) {
    return glm::log(obj);
}

/** Round an object of type T to the nearest integer. */
template <typename T>
auto round(T obj) {
    return glm::round(obj);
}

/** Multiply the elements of a vector together */
template <dimension D, typename T>
auto product(vec<D, T> vec) {
    auto result = (T)1;
    for (int d = 0; d < D; ++d) {
        result *= vec[d];
    }
    return result;
}

/** Sum the elements of a vector together */
template <dimension D, typename T>
auto sum(vec<D, T> vec) {
    auto result = (T)0;
    for (int d = 0; d < D; ++d) {
        result += vec[d];
    }
    return result;
}

/** Multiply the elements of a vector together */
template <dimension D>
int reduce(vec<D, int> vec) {
    int product = 1;
    for (int d = 0; d < D; ++d) {
        product *= vec[d];
    }
    return product;
}

/** Compute the norm of a vector of type T. */
template <dimension D, typename T>
T norm(vec<D, T> vec) {
    T squared_sum = 0;
    for (int d = 0; d < D; ++d) {
        squared_sum += vec[d] * vec[d];
    }
    return sqrt(squared_sum);
}

template <typename VecLike>
typename VecLike::value_type norm(VecLike x) {
    using T = typename VecLike::value_type;
    T result = (T)0;
    for (auto i = 0u; i < x.size(); ++i) {
        result += x[i] * x[i];
    }
    return sqrt(result);
}

template <typename VecLike>
typename VecLike::value_type dot(VecLike x, VecLike y) {
    using T = typename VecLike::value_type;
    T result = (T)0;
    for (auto i = 0u; i < x.size(); ++i) {
        result += x[i] * y[i];
    }
    return result;
}

/** Compute the square root of an object of type T. */
template <typename T>
constexpr auto sqrt(T obj) {
    return glm::sqrt(obj);
}

/** Normalize an object of type T. */
template <typename T>
auto normalize(T obj) {
    return glm::normalize(obj);
}

/** Return the value of the maximum component. */
template <dimension D, typename T>
T max_element(vec<D, T> v) {
    auto max = std::numeric_limits<T>::min();
    for (int d = 0; d < D; ++d) {
        if (v[d] > max) {
            max = v[d];
        }
    }
    return max;
}

/** Return the index of the maximum component. */
template <dimension D, typename T>
int max_index(vec<D, T> x) {
    int max_index = -1;
    auto max = std::numeric_limits<T>::min();
    for (int d = 0; d < D; ++d) {
        if (x[d] > max) {
            max = x[d];
            max_index = d;
        }
    }
    return max_index;
}

/** Compute the distance between two vectors. */
template <dimension D, typename T>
T distance(vec<D, T> a, vec<D, T> b) {
    return glm::distance(a, b);
}

/** Compute the distance between two generic objects. */
template <typename T>
T distance(T a, T b) {
    return abs(a - b);
}

/** Compute the inner-product of two vectors. */
template <dimension D, typename T>
T dot(vec<D, T> a, vec<D, T> b) {
    return glm::dot(a, b);
}

/** Compute the cross-product of two vectors. */
template <typename T>
vec3<T> cross(vec3<T> a, vec3<T> b) {
    return glm::cross(a, b);
}

/** Compute the 'cross-product' of two 2-dimensional vectors. */
template <typename T>
T cross(typename vec<2_D, T>::type a, typename vec<2_D, T>::type b) {
    return a.x * b.y - a.y * b.x;
}

template <tomo::dimension D, typename T>
auto array_to_vec(std::array<T, D> in) {
    tomo::math::vec<D, T> out;
    for (int d = 0; d < D; ++d) {
        out[d] = in[d];
    }
    return out;
}

template <tomo::dimension D, typename T>
auto vec_to_array(tomo::math::vec<D, T> in) {
    std::array<T, D> out;
    for (int d = 0; d < D; ++d) {
        out[d] = in[d];
    }
    return out;
}

/** Return the `axis`-th vector of the standard basis */
template <dimension D, typename T>
auto standard_basis(int axis) {
    auto result = vec<D, T>{};
    result[axis] = (T)1;
    return result;
}

/** Return a vector with the sign of each component as 1 or -1 */
template <dimension D, typename T>
vec<D, int> sign(vec<D, T> rhs) {
    auto result = vec<D, int>{};
    for (int d = 0; d < D; ++d) {
        if (approx_equal(rhs[d], (T)0)) {
            result[d] = 0;
        } else {
            result[d] = rhs[d] > 0 ? 1 : -1;
        }
    }
    return result;
}

} // namespace math
} // namespace tomo
