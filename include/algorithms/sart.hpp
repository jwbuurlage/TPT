#pragma once

#include <iostream>
#include <vector>

#include "../logging.hpp"

namespace tomo {

template <dimension D, typename T, class Geometry, class Projector>
image<D, T> sart(const volume<D>& v, const Geometry& g,
                 const sinogram<D, T, Geometry, Projector>& p,
                 double beta = 0.5, int iterations = 10,
                 optional<image<D, T>> initial = optional<image<D, T>>()) {
    image<D, T> f(v);
    if (initial)
        f = initial.value();
    Projector proj(v);

    // the size of a single block
    int k = g.groups()[0];

    // compute $w_i \cdot w_i$
    std::vector<T> w_norms(g.lines());
    int line_number = 0;
    for (auto line : g) {
        proj.reset(line);
        for (auto elem : proj) {
            w_norms[line_number] += elem.value * elem.value;
        }
        ++line_number;
    }

    auto f_next = f;

    for (int iter = 0; iter < iterations; ++iter) {
        int s = k;
        int t = 0;

        int row = 0;
        for (auto line : g) {
            proj.reset(line);
            if (s == k) {
                // we now update the image
                if (t > 0) {
                    f = f_next;
                }
                ++t;
                s = 0;
            }

            if (w_norms[row] > math::epsilon) {
                T alpha = 0.0;
                for (auto elem : proj) {
                    alpha += f[elem.index] * elem.value;
                }

                auto factor = beta * ((p[row] - alpha) / w_norms[row]);
                for (auto elem : proj)
                    f_next[elem.index] += factor * elem.value;
            }

            ++row;
            ++s;
        }
    }

    return f;
}

} // namespace tomo
