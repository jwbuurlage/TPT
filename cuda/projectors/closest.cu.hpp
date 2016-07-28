#pragma once

#include "../device_data.hpp"

namespace tomo {
namespace cuda {

// gpu, first just 'closest' kernel (value 1)
template <typename T, typename Func>
__device__ void project_closest(device::line<T> line, device::volume v, Func f) {
    device::vec2<T> current = line.origin;
    while (current.x <= v.x + EPSILON && current.y <= v.y + EPSILON &&
           current.x + EPSILON >= 0 && current.y + EPSILON >= 0) {
        int index = (int)(current.x) + v.x * (int)(current.y);
        if (index < v.x * v.y && index >= 0) {
            f(index);
        }

        current.x += line.delta.x;
        current.y += line.delta.y;
    }
}

} // namespace cuda
} // namespace tomo
