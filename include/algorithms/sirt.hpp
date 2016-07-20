#pragma once

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

    // first we compute R and C
    std::vector<T> R(g.lines());
    std::vector<T> C(v.cells());
    int line_number = 0;
    for (auto line : g) {
        Projector proj(line, v);
        for (auto elem : proj) {
            R[line_number] += elem.value;
            C[elem.index] += elem.value;
        }
        ++line_number;
    }

    for (auto& r : R)
        r = 1.0 / r;
    for (auto& c : C)
        c = 1.0 / c;

    sinogram<D, T, Geometry, Projector> s1(g);
    image<D, T> s2(v);

    for (int k = 0; k < iterations; ++k) {
        // next we compute Wx
        auto Wx = forward_projection(f, g);

        // compute R(p - Wx)
        for (int j = 0; j < g.lines(); ++j) {
            s1[j] = (p[j] - Wx[j]) * R[j];
        }
        // zero the image
        s2 = image<D, T>(v);

        // multiply with W^T
        int row = 0;
        for (auto line : g) {
            Projector proj(line, v);
            for (auto elem : proj) {
                s2[elem.index] += elem.value * s1[row];
            }
            ++row;
        }

        // update image while scaling with C
        for (int j = 0; j < v.cells(); ++j) {
            f[j] += beta * C[j] * s2[j];
        }
    }

    return f;
}

} // namespace tomo
