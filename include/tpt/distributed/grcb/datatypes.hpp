#pragma once

#include <array>

#include "tpt/tpt.hpp"

namespace tpt::grcb {

using namespace tpt::math;

template <typename T>
using cube = std::array<vec3<T>, 8>;

template <int D, typename T>
using box = std::array<vec<D, T>, 2>;

template <typename T>
using split_function =
    std::function<T(const geometry::base<3_D, T>&, cube<T>, size_t)>;

enum class flag {
    gradient_volume = 0x0,
    split_midway = 0x10,
    split_equal_shadow = 0x11,
    split_equal_load = 0x12,
    split_smart = 0x13,
    split_smart_simple = 0x14,
    precision_low = 0x20,
    precision_medium = 0x21,
    precision_high = 0x22,
    integration_two_dim = 0x30,
    integration_three_dim = 0x31,
    integration_cone = 0x32,
};

} // namespace tpt::grcb
