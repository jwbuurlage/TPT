#pragma once

#include <experimental/optional>
using std::experimental::optional;

#include "../common.hpp"
#include "../geometry.hpp"
#include "../image.hpp"
#include "../math.hpp"
#include "../sinogram.hpp"
#include "../volume.hpp"

namespace tomo {

template <dimension D, typename T, class Geometry, class Projector>
image<D, T> art(const volume<D>& v, const Geometry& g,
                const sinogram<D, T, Geometry, Projector>& p, double beta = 0.5,
                optional<image<D, T>> initial = optional<image<D, T>>()) {
    image<D, T> f(v);
    if (initial)
        f = initial.value();

    // compute $w_i \cdot w_i$
    std::vector<T> w_norms(g.lines());
    int line_number = 0;
    for (auto line : g) {
        Projector proj(line, v);
        for (auto elem : proj) {
            w_norms[line_number] += elem.value * elem.value;
        }
        ++line_number;
    }

    int m = 16;
    for (int k = 0; k < m; ++k) {
        int row = 0;
        for (auto line : g) {
            if (w_norms[row] < math::epsilon)
                continue;

            Projector proj(line, v);
            T alpha = 0.0;
            for (auto elem : proj)
                alpha += f[elem.index] * elem.value;

            auto factor = beta * ((p[row] - alpha) / w_norms[row]);
            for (auto elem : proj)
                f[elem.index] += factor * elem.value;

            ++row;
        }
    }

    return f;
}

} // namespace tomo
