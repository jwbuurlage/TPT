#pragma once

#include "constants.hpp"

#include <glm/glm.hpp>

namespace tomo {

using namespace tomo::literals;

namespace math {

/** Compute \f$a^n\f$. */
template <typename T>
constexpr T pow(T a, int n) {
    int result = a;
    // FIXME use exponential powering
    for (int i = 1; i < n; ++i) {
        result *= a;
    }
    return result;
}

template <typename T>
auto max(T a, T b) {
    return a > b ? a : b;
}

} // namespace math
} // namespace tomo
