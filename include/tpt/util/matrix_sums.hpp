#pragma once

#include "../common.hpp"
#include "../geometry.hpp"
#include "../math.hpp"
#include "../projector.hpp"

namespace tpt {

template <dimension D, typename T>
image<D, T> column_sums(const geometry::base<D, T>& geom,
                              dim::base<D, T>& kernel) {
    auto v = kernel.get_volume();
    auto result = image<D, T>(v);
    for (auto[idx, line] : geom) {
        (void)idx;
        for (auto elem : kernel(line)) {
            result[elem.index] += elem.value;
        }
    }

    return result;
}

template <dimension D, typename T>
projections<D, T> row_sums(const geometry::base<D, T>& geom,
                                 dim::base<D, T>& kernel) {
    auto result = projections<D, T>(geom);
    for (auto[idx, line] : geom) {
        for (auto elem : kernel(line)) {
            result[idx] += elem.value;
        }
    }

    return result;
}

} // namespace tpt
