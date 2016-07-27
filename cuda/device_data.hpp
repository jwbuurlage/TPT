#pragma once

namespace tomo {
namespace cuda {

template <typename T>
struct vec2 {
    T x;
    T y;

    void operator+=(vec2& other) {
        x += other.x;
        y += other.y;
    }
};

template <typename T>
struct device_line {
    vec2<T> origin;
    vec2<T> delta;
};

struct device_volume {
    int x;
    int y;
};

} // namespace cuda
} // namespace tomo
