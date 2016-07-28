#pragma once

namespace tomo {
namespace cuda {
namespace device {

#define EPSILON 1e-6

template <typename T>
struct vec2 {
    T x;
    T y;
};

template <typename T>
struct line {
    vec2<T> origin;
    vec2<T> delta;
};

struct volume {
    int x;
    int y;
};

} // namespace device
} // namespace cuda
} // namespace tomo
