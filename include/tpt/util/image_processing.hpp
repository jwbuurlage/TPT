#pragma once

#include "../image.hpp"
#include "../projections.hpp"
#include "../math.hpp"

namespace tpt {
namespace img {

template <dimension D, typename T>
image<D, T> clamp(image<D, T> img, T min, T max) {
    auto result = img;
    for (auto& x : result.mutable_data()) {
        if (x < min) {
            x = min;
        }
        if (x > max) {
            x = max;
        }
    }
    return result;
}

template <dimension D, typename T>
image<D, T> log(image<D, T> img) {
    auto result = img;
    for (auto& x : result.mutable_data()) {
        x = math::log(x);
    }
    return result;
}

template <dimension D, typename T>
image<D, T> operator-(image<D, T> img) {
    auto result = img;
    for (auto& x : result.mutable_data()) {
        x = -x;
    }
    return result;
}

template <dimension D, typename T>
image<D, T> operator+(const image<D, T>& lhs, const image<D, T>& rhs) {
    assert(lhs.size() == rhs.size());
    auto result = lhs;
    int idx = 0;
    for (auto& x : result.mutable_data()) {
        x += rhs[idx++];
    }
    return result;
}

template <dimension D, typename T>
image<D, T> operator*(const image<D, T>& lhs, const image<D, T>& rhs) {
    assert(lhs.size() == rhs.size());
    auto result = lhs;
    int idx = 0;
    for (auto& x : result.mutable_data()) {
        x *= rhs[idx++];
    }
    return result;
}

template <dimension D, typename T>
image<D, T> operator/(const image<D, T>& lhs, const image<D, T>& rhs) {
    assert(lhs.size() == rhs.size());
    auto result = lhs;
    int idx = 0;
    for (auto& x : result.mutable_data()) {
        x /= rhs[idx++];
    }
    return result;
}

template <dimension D, typename T>
image<D, T> operator*(T lhs, const image<D, T>& rhs) {
    auto result = rhs;
    for (auto& x : result.mutable_data()) {
        x *= lhs;
    }
    return result;
}

template <dimension D, typename T>
image<D, T>& operator+=(image<D, T>& lhs, const image<D, T>& rhs) {
    assert(lhs.size() == rhs.size());
    int idx = 0;
    for (auto& x : lhs.mutable_data()) {
        x += rhs[idx++];
    }
    return lhs;
}

template <dimension D, typename T>
projections<D, T> operator*(T lhs, const projections<D, T>& rhs) {
    auto result = rhs;
    for (auto& x : result.mutable_data()) {
        x *= lhs;
    }
    return result;
}

template <dimension D, typename T>
projections<D, T> operator-(const projections<D, T>& lhs, const projections<D, T>& rhs) {
    assert(lhs.size() == rhs.size());
    auto result = lhs;
    int idx = 0;
    for (auto& x : result.mutable_data()) {
        x -= rhs[idx++];
    }
    return result;
}

template <dimension D, typename T>
projections<D, T>& operator-=(projections<D, T>& lhs, const projections<D, T>& rhs) {
    assert(lhs.size() == rhs.size());
    int idx = 0;
    for (auto& x : lhs.mutable_data()) {
        x -= rhs[idx++];
    }
    return lhs;
}

template <dimension D, typename T>
image<D, T> operator-(const image<D, T>& lhs, const image<D, T>& rhs) {
    assert(lhs.size() == rhs.size());
    auto result = lhs;
    int idx = 0;
    for (auto& x : result.mutable_data()) {
        x -= rhs[idx++];
    }
    return result;
}

template <dimension D, typename T>
projections<D, T>& operator+=(projections<D, T>& lhs, const projections<D, T>& rhs) {
    assert(lhs.size() == rhs.size());
    int idx = 0;
    for (auto& x : lhs.mutable_data()) {
        x += rhs[idx++];
    }
    return lhs;
}

template <dimension D, typename T>
projections<D, T>& operator*=(projections<D, T>& lhs, T rhs) {
    for (auto& x : lhs.mutable_data()) {
        x *= rhs;
    }
    return lhs;
}

} // namespace img
} // namespace tpt
