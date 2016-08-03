#pragma once

#include "../device_data.hpp"

namespace tomo {
namespace cuda {

// gpu, first just 'closest' kernel (value 1)
template <typename T, typename Func>
__device__ void project_closest(device::line<T> line, device::volume v,
                                Func f) {
    device::vec2<T> current = line.origin;
    for (int i = 0; i < line.steps; ++i) {
        int index = (int)(current.x) + v.x * (int)(current.y);
        f(index);

        current.x += line.delta.x;
        current.y += line.delta.y;
    }
}

// gpu, first just 'closest' kernel (value 1)
template <typename T, typename Func>
__device__ void project_closest_texture(device::line<T> line, device::volume v,
                                        Func f) {
    device::vec2<T> current = line.origin;
    for (int i = 0; i < line.steps; ++i) {
        f(current.x, current.y);

        current.x += line.delta.x;
        current.y += line.delta.y;
    }
}

// gpu, first just 'closest' kernel (value 1)
template <typename T, typename Func>
__device__ void project_closest_surface(device::line<T> line, device::volume v,
                                        Func f) {
    device::vec2<T> current = line.origin;
    for (int i = 0; i < line.steps; ++i) {
        int x = (int)current.x;
        int y = (int)current.y;

        if (x < v.x && y < v.y && x >= 0 && y >= 0)
            f(x, y);

        current.x += line.delta.x;
        current.y += line.delta.y;
    }
}

} // namespace cuda
} // namespace tomo
