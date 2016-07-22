#pragma once

#include <experimental/optional>
using std::experimental::optional;

#include <vector>

namespace tomo {

template <dimension D, typename T, class Geometry, class Projector>
image<D, T> sirt(const volume<D>& v, const Geometry& g,
                 const sinogram<D, T, Geometry, Projector>& p,
                 double beta = 0.5, int iterations = 10,
                 optional<image<D, T>> initial = optional<image<D, T>>()) {
    image<D, T> f(v);
    if (initial)
        f = initial.value();
    Projector proj(v);

    // first we compute R and C
    std::vector<T> R(g.lines());
    std::vector<T> bC(v.cells());
    int line_number = 0;

    for (auto line : g) {
        proj.reset(line);
        for (auto elem : proj) {
            R[line_number] += elem.value;
            bC[elem.index] += elem.value;
        }
        ++line_number;
    }

    for (auto& r : R)
        r = 1.0 / r;
    for (auto& bc : bC)
        bc = beta / bc;

    sinogram<D, T, Geometry, Projector> s1(g);
    image<D, T> s2(v);
    for (int k = 0; k < iterations; ++k) {
        // compute Wx
        auto Wx = forward_projection(f, g, proj);

        // compute R(p - Wx)
        for (int j = 0; j < g.lines(); ++j) {
            s1[j] = (p[j] - Wx[j]) * R[j];
        }

        // zero the image
        s2.clear();

        // multiply with W^T
        int row = 0;
        for (auto line : g) {
            proj.reset(line);
            for (auto elem : proj) {
                s2[elem.index] += elem.value * s1[row];
            }
            ++row;
        }

        // update image while scaling with beta * C
        for (int j = 0; j < v.cells(); ++j) {
            f[j] += bC[j] * s2[j];
        }
    }

    return f;
}

} // namespace tomo
