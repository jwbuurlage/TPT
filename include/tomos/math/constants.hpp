#pragma once

#include <glm/gtc/constants.hpp>
#include <glm/glm.hpp>

namespace tomo {
namespace math {

/** A constant representing a small value of type T. */
template <typename T>
constexpr T epsilon = (T)1e-5;

/** A constant representing pi with the precision of type T. */
template <typename T>
constexpr auto pi = glm::pi<T>();

/** A constant representing sqrt(2) with the precision of type T. */
template <typename T>
constexpr T sqrt2 = (T)1.41421356237;

} // namespace math
} // namespace tomo
